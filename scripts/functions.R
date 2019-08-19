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
  , perm_path = NULL# path to save the permuted data sets and moving window correlations to
  , vhat_path = NULL# path to save the v hat values, if desired
  , out_path = NULL # path to save the output, if desired
  , ID = NULL       # option participant ID if saving data to an output file
){
  start_time <- Sys.time()
  # make sure required packages are loaded
  if(isNamespaceLoaded("parallel") == F){library(parallel)}
  if(isNamespaceLoaded("plyr") == F){library(plyr)}
  if(isNamespaceLoaded("tidyverse") == F){library(tidyverse)}
  # make sure data are sorted by time
  data <- unclass(data) %>% as.data.frame()
  data <- data[order(data[,timevar], na.last = NA, decreasing = F), ]
  
  # penalties for later
  v_max <- v_max_fun(data)
  penalty <- penalty_fun(kmax, v_max, nrow(data))
  
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
    filter(tau_p_min_1 != tau_p)%>% # make sure that the last obs the previous window is
    # before the current window 
    filter((tau_p_min_1 + 1) != tau_p &
           (tau_p_min_1 != (tau_p + 1))) # last obs of current phase must greater than last obs of previous phase
  
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
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores)
  
  # import global env variables for parallel computing
  clusterExport(cl, varlist = c("data", "v_hats_comb", "perm_all_fun", "mw_cor_fun", "v_max_fun", "v_hats",
                                "r_vec_fun", "v_hat_fun", "v_hat_par_fun", "gauss_kernal_fun", "jitter_fun",
                                "window", "start", "end", "kmax", "penalty_fun", "Mode"), envir=environment())
  clusterCall(cl, function() lapply(c("tidyverse", "psych"), require, character.only = TRUE))
  # run the variance function on raw data
  v_hats_comb$v_hat <- parallel::parLapply(cl, v_hats_comb$data, v_hat_par_fun)
  
  # run the variance function on permuted sets
  v_hats_comb_perm <- tibble(sample = 1:iter)
  v_hats_comb_perm$v_hat <- parLapply(cl, 1:iter,
                  function(x) perm_all_fun(data, window, start, end, kmax))
  # v_hats_comb_perm$v_hat <- lapply(1:iter, function(x) perm_all_fun(data, window, start, end, kmax))
  stopCluster(cl) # end parallel computing session
  
  # unnest the variances from the list
  v_hats_comb <- v_hats_comb %>% unnest(v_hat) %>% unnest(data)
  
  # join the variances for different phases back with all possible 
  # phases for different numbers of change points
  v_hats_comb <- v_hats_comb %>% #select(-rn) %>%
    full_join(v_hats) %>%
    arrange(k, comb)
  
  # extract the permuted data sets and moving window correlations
  # perm_data <- v_hats_comb_perm %>% 
  #   mutate(data = map(v_hat, ~.$data),
  #          mw_cor = map(v_hat, ~.$mw_cor)) %>%
  #   select(-v_hat)
  # if(!is.null(perm_path)){save(perm_data, file = sprintf("%s/perm_data_%s.Rdata", perm_path, ID))}
  # # get rid of permuted data sets to save memory
  # rm(perm_data)
  
  # pen_perm <- v_hats_comb_perm %>%
  #   mutate(penalty = map(v_hat, ~.$penalty)) %>%
  #   select(sample, penalty) %>%
  #   unnest(penalty)
  
  r_hats <- v_hats_comb %>%
    group_by(k, comb) %>%
    summarize(rhat = mean(v_hat, na.rm = T)) %>%
    ungroup()
  
  r_hats_perm <- v_hats_comb_perm %>%
    mutate(rhat = map(v_hat, ~.$r_hats)) %>%
    select(-v_hat)
  
  # v_hats_comb_perm <- v_hats_comb_perm %>% 
  #   mutate(v_hat = map(v_hat, ~.$v_hat))
  
  # if(!is.null(vhat_path)){save(v_hats_comb, v_hats_comb_perm, penalty, file = sprintf("%s/vhats_%s.Rdata", vhat_path, ID))}
  rm(v_hats_comb)
  # rm(v_hats_comb_perm)
  
  # extract the permuted variances and unnest them
  # then join the variances for different phases back with all possible 
  # phases for different numbers of change points
  # v_hats_comb_perm <- v_hats_comb_perm %>% 
  #   unnest(v_hat) %>% 
  #   unnest(v_hat) %>% 
  #   select(-data) %>%
  #   full_join(v_hats_comb %>% select(rn, tau_p, tau_p_min_1)) %>%
  #   select(-rn) %>%
  #   distinct() %>%
  #   mutate(scrap = 1) %>%
  #   full_join(v_hats %>% mutate(scrap = 1)) %>%
  #   select(-scrap) %>%
  #   arrange(k, comb, sample)
  
  # calculate the r_hat values (average within-phase variance)
  
  
  # calculate the r_hat values (average within-phase variance) for the permutations
  
  # r_hats_perm <- v_hats_comb_perm %>% 
  #   group_by(sample, k, comb) %>%
  #   summarize(rhat = mean(v_hat, na.rm = T)) %>%
  #   ungroup() 
  
  # rm(v_hats_comb_perm)
  
  # number of change points
  # variance test
  # P_{variancetest} = #(\hat{R}_{min,K=0,perm} > \hat{R}_{min, K=0}) / B
  # where B is the number of permutations
  
  r_hat_perm_k_0 <- r_hats_perm %>% 
    mutate(rhat = map(rhat, function(x) x %>% filter(k == 0))) %>%
    unnest(rhat)
  
  v_test <- sum(r_hat_perm_k_0$rhat > (r_hats %>% filter(k == 0))$rhat, na.rm = T) / iter
  
  if(v_test < .025){print("V_test: k > 0")} else {print("V_test: No change points")}
  
  # variance drop test
  # find minimized variances for each k
   r_min <- r_hats %>% 
    group_by(k) %>%
    mutate(mincomb = comb[rhat == min(rhat, na.rm = T)]) %>%
    filter(comb == mincomb) 
  
  # find drops between minimized variances for raw data
  v_drop_raw <- r_min %>% # minimized r_hats
    # summarize(raw_min = min(rhat, na.rm = T)) %>%
    ungroup() %>%
    mutate(raw_min_1 = lag(rhat), # minimized r_hats for k-1
           raw_drop = rhat - raw_min_1) %>% # difference in min rhat for k and k-1
    rename(raw_min = rhat)
  
  min_fun <- function(df){
    df %>%
      group_by(k) %>%
      mutate(mincomb = comb[rhat == min(rhat, na.rm = T)]) %>%
      filter(comb == mincomb)
  }
  
  r_min_perm <- r_hats_perm %>%
    mutate(r_min_perm = map(rhat, min_fun)) %>%
    unnest(r_min_perm, .drop = T)
  
  # # find minimized variances for each sample and k
  # r_min_perm <- r_hats_perm %>%
  #   filter(!is.na(sample)) %>%
  #   group_by(sample, k) %>%
  #   mutate(mincomb = comb[rhat == min(rhat, na.rm = T)]) %>%
  #   filter(comb == mincomb)
  
  # find drops between minimized variances for raw data
  v_drop_perm <- r_min_perm %>% # minimized r_hats
    group_by(sample) %>%
    mutate(min_1 = lag(rhat), # minimized r_hats for k-1
           drop = rhat - min_1) %>% # difference in min rhat for k and k-1
    ungroup() %>%
    rename(min = rhat)
  
  # join raw and permuted variance drops
  v_drop <- v_drop_perm %>%
    select(sample, k, min, min_1, drop) %>%
    full_join(v_drop_raw %>% select(-comb, -mincomb)) 
  
  # now for each combo of k, match the permutated with the raw drops
  v_drop_test <- v_drop %>%
    filter(k != 0) %>%
    group_by(sample) %>%
    summarize(max_drop = max(abs(drop), na.rm = T), # find the max drop across all k's and samples
              max_raw_drop = max(abs(raw_drop), na.rm = T)) %>% # find the max raw drop across all k's
    ungroup() %>%
    # proportion of max permuted drops that are greater than raw drops 
    # (they shouldn't be if there are real changes)
    summarize(p = sum(max_drop > max_raw_drop, na.rm = T)/iter) 
  
  # print results of the Variance drop test
  if(any(v_drop_test$p < .025)){print("V_drop: Change Points detected")} else {print("V_drop: No Change Points")}
  
  
  k <- 0 # empirical k, will be changed if either tests are significant
  knots <- NA # empirical location of change points, will be changed if k > 0
  C <- 1 # penalization, will be selected is k > 0
  if(v_drop_test$p < .025 | v_test < .025){
    # choosing k, from Cabrieto et al (2018b), Information Sciences paper
    k_raw <- r_hats %>% # penalize rhat based on k and a constant C
      # filter(k != 0) %>%
      # bring in the penalties for all rhats
      full_join(penalty %>% full_join(crossing(k=0, C=1:500, penalty = 0))) %>%
      # if k == 0, there is no penalty
      mutate(penalty = ifelse(is.na(penalty), 0, penalty),
             raw_khat = rhat + penalty) # penalized averge within-phase variance
    
    # choose the penalty using grid search
    k_choose <- k_raw %>% group_by(C) %>%
      # find the combination with the minimized khat value for each penalization (C)
      # as well as the khat value for that combination
      mutate(min_comb = comb[raw_khat == min(raw_khat, na.rm = T)],
             min_k = k[raw_khat == min(raw_khat, na.rm = T)]) %>%
      ungroup() %>%
      # keep only the minimized khat values
      filter(comb == min_comb & k == min_k) %>%
      # get rid of zeros because they are just here to help us choose C
      # C = first C where k == 0 is chosen - 1
      filter(k != 0)
  
    mode_k <- Mode(k_choose$k) # choose k with the most stable mode for k
    comb_k <- (k_choose %>% filter(k == mode_k) %>% arrange(C))$comb # save the combination that matches
    
    C <- min((k_choose %>% filter(k == mode_k))$C) - 1 # penalization
    k <- mode_k
    knots <- (v_hats %>% filter(k == mode_k & comb == comb_k))$tau_p_min_1 + 1
  }
  
  # save change point results in a list
  cp <- list(
    cp = ifelse((v_test < .025 | v_drop$p < .025), "CP", "NCP"),
    knots = knots,
    k = k,
    C = C
    )
  
  # create list of all results
  out = list(
    v_drop = v_drop
    , v_test = v_test
    , r_hats = r_hats
    , r_hats_perm = r_hats_perm
    , cp = cp
  )
  
  # save the results to file or return as an object
  if(!is.null(out_path)){
    save(out, file = sprintf("%s/results_%s.Rdata", out_path, ID))
    return(T)
  } else {
      return(out)
  }
  end_time <- Sys.time()
  print(paste("Time elapsed", round(end_time-start_time,2)))
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

perm_all_fun <- function(data, window, start, end, k_max){
  n <- nrow(data) # number of observations
  s_data <- sample_n(data, n, replace = F) # shuffle the data
  v_max <- v_max_fun(s_data)
  # penalty <- penalty_fun(k_max, v_max, nrow(s_data))
  
  mw_cor <- mw_cor_fun(s_data, window = window) # moving window correlations
  mw_r_nested <- mw_cor %>% select(i, r) %>% group_by(i) %>% nest()
  
  # crossed combinations of starting and ending points
  v_hats_comb <- crossing(
    tau_p = start:end,
    tau_p_min_1 = start:end
  ) %>%
    filter(tau_p_min_1 != tau_p) %>%# make sure that the last obs the previous window is
  # before the current window 
    filter((tau_p_min_1 + 1) != tau_p &
           (tau_p_min_1 != (tau_p + 1)))
  
  # match the starting observations with the correct windows of correlations
  v_hats_comb <- v_hats_comb %>%
    mutate(rn = 1:n(), 
           R = list(mw_r_nested)) %>%
    group_by(rn) %>%
    nest()
  
  # get the variance of the moving window correlations within a window
  v_hats_comb$v_hat <- lapply(v_hats_comb$data, v_hat_par_fun)
  
  v_hats_comb <- v_hats_comb %>% unnest(v_hat) %>% unnest(data)
  
  # join the variances for different phases back with all possible 
  # phases for different numbers of change points
  v_hats_comb <- v_hats_comb %>% #select(-rn) %>%
    full_join(v_hats) %>%
    arrange(k, comb)
  
  r_hats <- v_hats_comb %>%
    group_by(k, comb) %>%
    summarize(rhat = mean(v_hat, na.rm = T)) %>%
    ungroup()
  
  # return the results
  return(list(
    r_hats = r_hats
    # , data = s_data
    # , mw_cor = mw_cor
    # , v_hat = v_hats_comb
    # , v_max = v_max
  )
    )
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
    filter(inner_l != outer_l) %>%
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

v_max_fun <- function(
  data
) {
  # first 5% of obs
  d_beg_5 <- data %>% filter(week %in% 1:ceiling(quantile(1:nrow(data), .05))) %>% select(-week)
  # last 5% of obs
  d_las_5 <- data %>% filter(week %in% (floor(quantile(1:nrow(data), .95))):nrow(data)) %>% select(-week)
  # calculate the trace of the covariance matrices of first and last 5%
  v_beg <- tr(cov(d_beg_5))
  v_las <- tr(cov(d_las_5))
  v_max <- ifelse(v_beg > v_las, v_beg, v_las)
}

penalty_fun <- function(
  k_max
  , v_max
  , n
) {
  pen_math_fun <- function(C, k, n, v_max){
    pen_k <- C*((v_max*(k + 1))/n)*(1 + log((n/(k_max+1))))
  }
  
  penalties <- crossing(
    C = 1:500,
    k = 1:k_max,
    n = n,
    v_max = v_max
  ) %>%
    mutate(penalty = pmap_dbl(list(C, k, n, v_max), pen_math_fun))
}

Mode <- function(x) {
  ux <- unique(x)
  ux <- ux[!is.na(ux)]
  ux[which.max(tabulate(match(x, ux)))]
}

vhats_setup_fun <- function(kmax, start, end){
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
  return(v_hats)
}

load_fun <- function(sid){
  load(sprintf("%s/out/results_%s.RData", "~/Downloads", sid))
  data <- (KCP_nested %>% filter(ID == sid))$data[[1]] 
  start <- (1 + 10 %/% 2) # e.g. if the window was 10, the start will be 5
  end <- (nrow(data) - 10 %/% 2) # e.g. window = 10 and nrow = 56, end will be 51
  v_max <- v_max_fun(data)
  penalty <- penalty_fun(3, v_max, nrow(data))
  v_hats <- vhats_setup_fun(3, start, end)
  k <- 0 # empirical k, will be changed if either tests are significant
  knots <- NA # empirical location of change points, will be changed if k > 0
  r_hats <- out$r_hats
  iter <- 500
  v_drop_test <- out$v_drop %>%
    filter(k != 0) %>%
    group_by(sample) %>%
    summarize(max_drop = max(abs(drop), na.rm = T), # find the max drop across all k's and samples
              max_raw_drop = max(abs(raw_drop), na.rm = T)) %>% # find the max raw drop across all k's
    ungroup() %>%
    # proportion of max permuted drops that are greater than raw drops 
    # (they shouldn't be if there are real changes)
    summarize(p = sum(max_drop > max_raw_drop, na.rm = T)/iter) 
  v_test <- out$v_test
  out$v_drop_test <- v_drop_test
  C <- 1 # penalization, will be selected is k > 0
  if(v_drop_test$p < .025 | v_test < .025){
    # choosing k, from Cabrieto et al (2018b), Information Sciences paper
    k_raw <- r_hats %>% # penalize rhat based on k and a constant C
      # filter(k != 0) %>%
      # bring in the penalties for all rhats
      full_join(penalty %>% full_join(crossing(k=0, C=1:500, penalty = 0))) %>%
      # if k == 0, there is no penalty
      mutate(penalty = ifelse(is.na(penalty), 0, penalty),
             raw_khat = rhat + penalty) # penalized averge within-phase variance
    
    # choose the penalty using grid search
    
    k_choose <- k_raw %>% 
      left_join(v_hats %>% filter(tau_p == tau_p_min_1 | (tau_p == end & tau_p_min_1 %in% (end-1:4))) %>% select(k, comb) %>% mutate(bad = "yes")) %>% 
      group_by(C) %>%
      # find the combination with the minimized khat value for each penalization (C)
      # as well as the khat value for that combination
      mutate(min_comb = comb[raw_khat == min(raw_khat, na.rm = T)],
             min_k = k[raw_khat == min(raw_khat, na.rm = T)]) %>%
      ungroup() %>%
      # keep only the minimized khat values
      filter(comb == min_comb & k == min_k) %>%
      # get rid of zeros because they are just here to help us choose C
      # C = first C where k == 0 is chosen - 1
      filter(k != 0)
    
    mode_k <- Mode(k_choose$k) # choose k with the most stable mode for k
    comb_k <- (k_choose %>% filter(k == mode_k) %>% arrange(C))$comb # save the combination that matches
    
    C <- min((k_choose %>% filter(k == mode_k))$C) - 1 # penalization
    k <- mode_k
    knots <- (v_hats %>% filter(k == mode_k & comb == comb_k & tau != "i"))$tau_p_min_1 -1
  }
  out$cp$knots <- knots; out$cp$k <- k; out$cp$C <- C
  out <- out[which(names(out) %in% c("cp", "v_drop", "v_test"))]
  return(out)
}