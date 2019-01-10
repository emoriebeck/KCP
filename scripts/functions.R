# library(astrochron)
# library(parallel)
# library(plyr)
# library(tidyverse)

# ex <- cycles(freqs=c(1/40,1/20),noisevar=.2)
# ex[3] <- cycles(freqs=c(1/40,1/20),noisevar=0.2)[2]
# ex[1]=ex[1]+rnorm(500,sd=5)
# ex = ex[order(ex[1],na.last=NA,decreasing=FALSE),]
# data <- ex[1:65,]
# kmax <- 2
# window <- 5
# timevar <- "time"
# iter <- 100


KCP <- function(
  data              # N x P numeric matrix or data frame
  , kmax = 1        # max number of change points to test
  , window = 5      # the size of the moving window used for correlations
  , timevar         # name of the time variable for sorting
  , iter = 1000     # number of permutations
  , alpha = .05     # "significance criterion"
  , rule = "OR"     # should both variance and variance drop tests be required to be significant
  , perm_data = NULL# path to save the permuted data sets and moving window correlations to
  , vhat_data = NULL# path to save the v hat values, if desired
  , ID = NULL       # option participant ID if saving data to an output file
){
  # make sure required packages are loaded
  if(isNamespaceLoaded("parallel") == F){library(parallel)}
  if(isNamespaceLoaded("plyr") == F){library(plyr)}
  if(isNamespaceLoaded("tidyverse") == F){library(tidyverse)}
  # make sure data are sorted by time
  data <- unclass(data) %>% as.data.frame()
  data <- data[order(data[,timevar], na.last = NA, decreasing = F), ]
  start <- (1 + window %/% 2) # e.g. if the window was 10, the start will be 5
  end <- (nrow(data) - window %/% 2) # e.g. window = 10 and nrow = 56, end will be 51
  
  # calculate the moving window correlations for the raw data
  raw_mw_r <- mw_cor_fun(data %>% select(-one_of(timevar)), window = window)
  raw_mw_r_nested <- raw_mw_r %>% select(i, r) %>% group_by(i) %>% nest()
  
  # calculate all combinations of possible change points in the time series
  # there will be overlap across different sets, so we are just setting this 
  # up for later for now. Later, we'll join it with the calculated within-phase
  # variances for further analysis
  if(kmax >= 0){
    # no change points
    v_hats <- crossing(
      k = 0, # k = 0 means no change points
      tau = "i", # i = initial time series point
      tau_p_min_1 = start, # beginning of previous phase
      tau_p = end, # last obs of current phase
      comb = 1 # needed later for keeping things together when there are change points
    ) 
  } 
  if(kmax >= 1){
    v_hats <- crossing(
      k = 1, # k = 1 means 1 change point
      i = start, # setting last phase of initial phase to the median of the first time point
      a = (start + 1):end # set change point to all possible locations between the beginning and end of the series
    ) %>%
      filter(a > i) %>% # last obs of current phase must be larger than last obs of last phase
      mutate(comb = 1:n()) %>% # needed later for keeping things together when there are change points
      gather(key = tau, value = v, -k, -comb) %>% # change to long
      group_by(comb) %>% # group by phases
      mutate(lead = lead(v)) %>% # shift last obs up one row to match with previous phases
      ungroup() %>%
      mutate(tau_p = ifelse(tau == "a", end, lead - 1), # tau_p is the last obs of current phase
             tau_p_min_1 = ifelse(tau == "i", start, v)) %>% # tau_p_min_1 is the last obs of the prvious phase
      select(k, comb, tau, tau_p, tau_p_min_1) %>% # keep only necessary columns
      full_join(v_hats) # join with previous v_hat combinations for fewer phases
  } 
  if(kmax >= 2){
    v_hats <-
          crossing(
            k = 2,  # k = 2 means 2 change points
            i = start,# setting last phase of initial phase to the median of the first time point
            a = (start + 1):(end), # set change point to all possible locations between the beginning and end of the series
            b = (a + 1):(end) # set change point to all possible locations between a and end of the series
          ) %>%
        filter((a > i) & (b > a)) %>% # last obs of current phase must be larger than last obs of last phase
        mutate(comb = 1:n()) %>% # needed later for keeping things together when there are change points
        gather(key = tau, value = v, -k, -comb) %>% # change to long
        group_by(comb) %>% # group by phases
        mutate(lead = lead(v)) %>% # shift last obs up one row to match with previous phases
        ungroup() %>%
        mutate(tau_p = ifelse(tau == "b", end, lead - 1), # tau_p is the last obs of current phase
               tau_p_min_1 = ifelse(tau == "i", start, v)) %>% # tau_p_min_1 is the last obs of the prvious phase
        select(k, comb, tau, tau_p, tau_p_min_1) %>% # keep only necessary columns
        full_join(v_hats) # join with previous v_hat combinations for fewer phases
    } 
  if (kmax == 3) {
      v_hats <- 
        crossing(
          k = 3, # k = 3 means 3 change points
          i = start,# setting last phase of initial phase to the median of the first time point
          a = (start + 1):(end),  # set change point to all possible locations between the beginning and end of the series
          b = (a + 1):(end), # set change point to all possible locations between a and end of the series
          c = (b + 1):end  # set change point to all possible locations between b and end of the series
        ) %>%
        filter((a > i) & (b > a) & (c > b)) %>% # last obs of current phase must be larger than last obs of last phase
        mutate(comb = 1:n()) %>% # needed later for keeping things together when there are change points
        gather(key = tau, value = v, -k, -comb) %>% # change to long
        group_by(comb) %>% # group by phases
        mutate(lead = lead(v)) %>% # shift last obs up one row to match with previous phases
        ungroup() %>%
        mutate(tau_p = ifelse(tau == "c", end, lead - 1), # tau_p is the last obs of current phase
               tau_p_min_1 = ifelse(tau == "i", start, v)) %>% # tau_p_min_1 is the last obs of the prvious phase
        select(k, comb, tau, tau_p, tau_p_min_1) %>% # keep only necessary columns
        full_join(v_hats) # join with previous v_hat combinations for fewer phases
    } 
  if (kmax > 3) {
      stop ("these are psychological time series and > 3 CP's isn't likely plausible")
  }
  
  # real data
  v_hats_comb <- crossing(
    tau_p = start:end, # last obs of current phase
    tau_p_min_1 = start:end # last obs of previous phase
  ) %>%
    filter((tau_p_min_1+2) < tau_p) # last obs of current phase must greater than last obs of previous phase
  
  # v_hats_perm <- v_hats %>% mutate(scrap = 1) %>% 
  #   arrange(k, comb) %>%
  #   full_join(tibble(sample = 1:iter, scrap = 1))
  
  v_hats_comb <- v_hats_comb %>%
    mutate(rn = 1:n(), # saving this row number for matching stuff later
           R = list(raw_mw_r %>% select(i,r) %>% group_by(i) %>% nest())) %>% # add in the moving window correlations
    group_by(rn) %>%
    nest() # nest for later use
  
  
  # function for running the permutations in parallel
  # given a list that contains tau_p, tau_p_min_1 and moving correlations, 
  # call the variance function to calculate variance in similarity
  v_hat_par_fun <- function(q){
    v_hat_fun(q$tau_p, q$tau_p_min_1, q$R[[1]])
  }
  
  # Calculate the number of cores
  no_cores <- detectCores() - 2
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  # import global env variables for parallel computing
  clusterExport(cl, varlist = c("data", "v_hats_comb", "perm_all_fun", "mw_cor_fun",
      "r_vec_fun", "v_hat_fun", "v_hat_par_fun", "gauss_kernal_fun", "jitter_fun",
      "window", "start", "end"))
  clusterCall(cl, function() library(tidyverse))
  # run the variance function on raw data
  v_hats_comb$v_hat <- parallel::parLapply(cl, v_hats_comb$data, v_hat_par_fun)
  
  # run the variance function on permuted sets
  v_hats_comb_perm <- tibble(sample = 1:iter)
  v_hats_comb_perm$v_hat <- parLapply(cl, 1:iter,
          function(x) perm_all_fun(data, window, start, end))
  # v_hats_comb_perm$v_hat <- lapply(1:iter,
  #           function(x) perm_all_fun(data, window, start, end))
  
  stopCluster(cl) # end parallel computing session
  
  # unnest the variances from the list
  v_hats_comb <- v_hats_comb %>% unnest(v_hat) %>% unnest(data)
  
  # join the variances for different phases back with all possible 
  # phases for different numbers of change points
  v_hats_comb <- v_hats_comb %>% #select(-rn) %>%
    full_join(v_hats) %>%
    arrange(k, comb)
  
  # extract the permuted data sets and moving window correlations
  perm_data <- v_hats_comb_perm %>% 
    mutate(data = map(v_hat, ~.$data),
           mw_cor = map(v_hat, ~.$mw_cor)) %>%
    select(-v_hat)
  if(!is.null(perm_data)){save(perm_data, file = sprintf("%s/perm_data_%s.Rdata", perm_data, ID))}
  
  # get rid of permuted data sets to save memory
  v_hats_comb_perm <- v_hats_comb_perm %>% 
    mutate(v_hat = map(v_hat, ~.$v_hat))
  rm(perm_data)
  
  # extract the permuted variances and unnest them
  # then join the variances for different phases back with all possible 
  # phases for different numbers of change points
  v_hats_comb_perm <- v_hats_comb_perm %>% 
    unnest(v_hat) %>% 
    unnest(v_hat) %>% 
    select(-data) %>%
    full_join(v_hats_comb %>% select(rn, tau_p, tau_p_min_1)) %>%
    select(-rn) %>%
    mutate(scrap = 1) %>%
    full_join(v_hats %>% mutate(scrap = 1)) %>%
    select(-scrap)
  
  if(!is.null(vhat_data)){save(v_hats_comb, v_hats_comb_perm, file = sprintf("%s/vhats_%s.Rdata", vhat_data, ID))}
  
  # calculate the r_hat values (average within-phase variance)
  r_hats <- v_hats_comb %>%
    group_by(k, comb) %>%
    summarize(rhat = mean(v_hat, na.rm = T)) %>%
    ungroup()
  
  # calculate the r_hat values (average within-phase variance) for the permutations
  r_hats_perm <- v_hats_comb_perm %>% 
    group_by(sample, k, comb) %>%
    summarize(rhat = mean(v_hat, na.rm = T)) %>%
    ungroup()
  
  # number of change points
  # variance test
  # P_{variancetest} = #(\hat{R}_{min,K=0,perm} > \hat{R}_{min, K=0}) / B
  # where B is the number of permutations
  v_test <- sum((r_hats_perm %>% filter(k == 0))$rhat > (r_hats %>% filter(k == 0))$rhat, na.rm = T) / iter
  
  if(v_test < .05){print("k > 0")} else {print("No change points")}
  
  # variance drop test
  v_drop_raw <- r_hats %>% 
    group_by(k, comb) %>%
    summarize(raw_min = min(rhat, na.rm = T))
  
  v_drop_perm <- r_hats_perm %>%
    filter(!is.na(sample)) %>%
    group_by(sample, k, comb) %>%
    summarize(min = min(rhat, na.rm = T)) %>%
    ungroup() %>%
    full_join(v_drop_raw)
  
  # now for each combo of k, match the permutated with the raw drops
  v_drop <- v_drop_perm %>%
    group_by(k, comb) %>%
    summarize(total = sum(min > raw_min, na.rm = T))
  
  out = list(
    v_drop = v_drop
    , v_test = v_test
    , r_hats = r_hats
    , r_hats_perm = r_hats_perm
  )
}

### mw_cor_fun ###

# function meant to take in data and call r_vec_fun to get moving window
# correlations across a set window size

mw_cor_fun <- function(
  data # shuffled or raw data matrix
  , window = 5 # the size of the moving window used for correlations
) {
  # setup
  ipts <- nrow(data) # number of obs
  nwin <- ipts - (window %/% 2) # get number of windows: N - window size + 1
  
  # get vectors of pairwise correlations across moving windows
  r_vec <- tibble(i = (1 + window %/% 2):nwin) %>% # data frame of just 1 to number of windows
    mutate(r = map(i, ~r_vec_fun(data, ., window))) %>% # create nested mw correlations
    unnest(r) # unnest to long format
}


### r_vec_fun ###

# Function used to calculate pairwise correlations for a single window 
# and save them in long format for later use as vectors

r_vec_fun <- function(
  data # N x P data matrix of raw or shuffled data
  , median # starting point of window
  , window # window size
  ) {
  # setup 
  start <- median - (window %/% 2)
  data <- data[start:(start+window-1),] # keep only data in window
  data <- jitter_fun(data)
  # run and clean correlations
  r <- cor(data, use = "pairwise") # calculate pairwise correlations
  r[upper.tri(r, diag = T)] <- NA # remove diagonal and duplicate correlations
  r %>% data.frame %>% # convert to a data frame
    mutate(V1 = rownames(.)) %>% # save the variable names
    gather(key = V2, value = r, -V1, na.rm = T) # change to long format
}

### perm_all_fun ###

# function meant to take in sorted data matrix, pull permutation samples, 
# and run mw correlations for each permutation
# returns a nested data frame

perm_all_fun <- function(data, window, start, end){
  n <- nrow(data) # number of observations
  s_data <- sample_n(data, n, replace = F) # shuffle the data
  mw_cor <- mw_cor_fun(s_data, window = window) # moving window correlations
  mw_r_nested <- mw_cor %>% select(i, r) %>% group_by(i) %>% nest()
  
  # crossed combinations of starting and ending points
  v_hats_comb <- crossing(
    tau_p = start:end,
    tau_p_min_1 = start:end
  ) %>%
    filter((tau_p_min_1+2) < tau_p) # make sure that the last obs the previous window is
  # before the current window
  
  # match the starting observations with the correct windows of correlations
  v_hats_comb <- v_hats_comb %>%
    mutate(rn = 1:n(), 
           R = list(mw_r_nested)) %>%
    group_by(rn) %>%
    nest()
  
  # get the variance of the moving window correlations within a window
  v_hats_comb$v_hat <- lapply(v_hats_comb$data, v_hat_par_fun)
  # return the results
  return(list(data = s_data, mw_cor = mw_cor, v_hat = v_hats_comb))
}

gauss_kernal_fun <- function(
  R_i # vector of correlations at i
  , R_j # vector of correlations at j
) {
  
  # calculate h
  # h is the median Euclidean distance among all X_i's
  h <- median(dist(R_i,method="euclidean"))
  
  # Gaussian similarity
  # Gk(R_i, R_j) = exp( (-||R_i - R_j||^2) / 2*h_R^2)
  num <- -1 * (norm(as.matrix(R_i - R_j))^2)
  den <- 2 * h^2
  Gk <- exp(num/den)
}

v_hat_fun <- function(
  tau_p # last obs in current phase
  , tau_p_min_1 # last obs in previous phase
  # , data # N x P data matrix of raw or shuffled data
  , R # long data frame of mw correlations with different median phases, i
) {
  start <- tau_p_min_1 # starting value
  # data <- data[start:tau_p,] # filter the data?
  outer <- (1 / (tau_p - tau_p_min_1)) # outer part of problem
  
  print(paste(tau_p, tau_p_min_1, sep = " "))
  inner <- crossing(inner_l = (tau_p_min_1+1):tau_p,
           outer_l = (tau_p_min_1+1):tau_p) %>%
    filter(inner_l != outer_l & inner_l < outer_l) %>%
    left_join(R %>% dplyr::select(inner_l = i, R_i = data)) %>%
    left_join(R %>% dplyr::select(outer_l = i, R_j = data)) %>%
    mutate(inner = map2_dbl(R_i, R_j, possibly(gauss_kernal_fun, NA_real_)))
 
  inner_sum <- sum(inner$inner, na.rm = T) # sum inner loop
  res <- (tau_p - tau_p_min_1) - (outer * inner_sum) # multiply outer by inner loop sum
}

jitter_fun <- function(df){
  sd_fun <- function(x){if(sd(x, na.rm = T) == 0) jitter(x, amount = runif(1,0,.05)) else x}
  df2 <- data.frame(apply(df, 2, sd_fun))
  colnames(df2) <- colnames(df2)
  return(df2)
}
