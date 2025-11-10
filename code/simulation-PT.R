library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(sandwich)
library(lmtest)
library(broom)
library(osqp)
library(Matrix)
library(CVXR)
library(synthdid)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("all-func.R")
set.seed(0)

# read repeated cross-sectional data; see 'code/prep-data.R' for how this data is generated
rcs_CPS <- read.csv("cps_rep_cross_sec.csv")
data(CPS) # load the CPS data provided by the synthdid package
CPS$state_num <- as.numeric(CPS$state)

# get (kt)-cell means and sd
cell_stats <- rcs_CPS %>%
  group_by(state_num, year) %>%
  summarise(
    n_kt=n(),
    wage=mean(outcome, na.rm = TRUE),
    sd=sd(outcome, na.rm = TRUE),
    .groups = "drop"
  )

source("synthdid-sdid-paper/R/placebo-simulations.R") # functions used to construct simulation DGPs in AAHIW
# mean outcome: wage, in N by T matrix form
cell_stats <- left_join(cell_stats, CPS[, c("state_num", "year", "min_wage")], by=c("state_num", "year"))
Y.wage = panel.matrices(data.frame(cell_stats), treatment='min_wage', outcome='wage', treated.last=FALSE)$Y

# assignment variable: whether minimum wage is higher than federal minimum wage in year 2000
last.col = function(X) { X[, ncol(X)] }
w.minwage = last.col(panel.matrices(data.frame(cell_stats), treatment='min_wage', treated.last=FALSE)$W)

# the following function from synthdid-sdid-paper/R/placebo-simulations.R is modified to return the assignment information
simulate_dgp = function(parameters, N1, T1){
  F_=parameters$F
  M=parameters$M
  Sigma = parameters$Sigma
  pi = parameters$pi
  
  N <- nrow(M)
  T_ <- ncol(M)
  
  assignment <- randomize_treatment(pi,N,N1)
  N1 <- sum(assignment)
  N0 <- N - N1
  T0 <- T_ - T1
  
  L <- F_ + M
  E <- mvtnorm::rmvnorm(N, sigma = Sigma)
  Y =  L + E 
  #### ------ START OF MODIFICATION ------
  return(list(Y=Y[order(assignment),], L=L[order(assignment),], E=E[order(assignment),], N0=N0, T0=T0, assignment=assignment))
  #### ------ END OF MODIFICATION ------ 
  # order units by treatment, so treated units are at the bottom of our data matrix in accordance with our convention
}

# fit a rank-4 factor matrix L to the outcome; decompose L = F + M; estimate error AR2 covariance E; estimate assignment probability
factor_model=estimate_dgp(Y.wage, w.minwage, rank=4)

# simulation configuration
true_tau <- 0 # null effect
num_MC <- 500
cand_tau_size <- 500
num_bstp_rep <- 1000
cand_tau_PT <- NULL
cand_tau_noPT <- NULL
frac_rej_PT <- c()
frac_rej_noPT <- c()
violation_noPT <- c()

# collect results; type in c("est", "se", "coverage", "cs_lb", "cs_ub", "time", "est_lb", "est_ub"); DGP in c(1,2)
result_pt <- data.frame(matrix(ncol=42, nrow = 0))
colnames(result_pt) <- c(1979:2016, 2018, "type", "DGP", "niter")
result_pt <- result_pt %>%
  mutate(type = as.character(type))

result_sdid <- data.frame(matrix(ncol=4, nrow = 0))
colnames(result_sdid) <- c(2018, "type", "DGP", "niter")
result_sdid <- result_sdid %>%
  mutate(type = as.character(type))

result_spt <- data.frame(matrix(ncol=4, nrow = 0))
colnames(result_spt) <- c(2018, "type", "DGP", "niter")
result_spt <- result_spt %>%
  mutate(type = as.character(type))

result_sc <- data.frame(matrix(ncol=4, nrow = 0))
colnames(result_sc) <- c(2018, "type", "DGP", "niter")
result_sc <- result_sc %>%
  mutate(type = as.character(type))


#### ---- simulation starts -----
for (iter in 1:num_MC){
  print(paste0("iter:  ", iter))
  start <- Sys.time()
  ### START: DGP-1 treated unit is one state, PT is satisfied ----
  # simulate Y = L + true_tau*W + E for SDID
  factor_model_1tr <- simulate_dgp(factor_model, N1=1, T1=1)
  
  # get the treated unit index
  trt_idx <- factor_model_1tr$assignment
  
  # generate treated unit's mean outcomes, trends, and control weights under PT
  trt_unit <- gen_trt_unit(trt_idx=trt_idx)
  pi_0 <- trt_unit$pi_0
  pi_1 <- trt_unit$pi_1
  
  # generate micro-level
  n_1 <- sum(subset(cell_stats, year==2018 & state_num %in% which(trt_idx==1))$n_kt) # total number of treated individuals
  n_0 <- sum(subset(cell_stats, year==2018 & state_num %in% which(trt_idx==0))$n_kt) # total number of control individuals
  n_k0 <- ceiling(n_0*pi_0) # control individuals sampled proportional to pi_0
  n_k1 <- ceiling(n_1*pi_1) # treated individuals sampled proportional to pi_1
  
  # repeated cross-sections for the control units
  rcs_data_ctrl <- map2_dfr(which(trt_idx==0), n_k0, \(k, nk)
                            map_dfr(1979:2018, \(t)
                                    gen_micro_data_control(k, t, nk)))
  
  # repeated cross-sections for the treated unit
  rcs_data_trt <- map2_dfr(which(trt_idx==1), n_k1, \(k, nk)
                           map_dfr(1979:2018, \(t)
                                   gen_micro_data_trt(trt_info=trt_unit,
                                                      k=k, t=t, n=nk,
                                                      true_tau=true_tau)))
  
  # combine into one data frame
  rcs_data <- bind_rows(rcs_data_ctrl, rcs_data_trt)
  
  # set T_0=2017 to be the reference period
  rcs_data$T_f <- factor(rcs_data$T_i, 
                         levels = sort(unique(rcs_data$T_i)))
  rcs_data$T_f <- relevel(rcs_data$T_f, ref = "2017")  
  
  ## START: PT-DID ------------------------------------
  # double check 2018 coeff matches difference-in-means estimator
  trt_2018 <- subset(rcs_data, T_i==2018 & D_i==1)$Y_i
  trt_2017 <- subset(rcs_data, T_i==2017 & D_i==1)$Y_i
  ctrl_2018 <- subset(rcs_data, T_i==2018 & D_i==0)$Y_i
  ctrl_2017 <- subset(rcs_data, T_i==2017 & D_i==0)$Y_i
  diff_in_means <- mean(trt_2018)-mean(trt_2017)-(mean(ctrl_2018)-mean(ctrl_2017))
  
  # difference-in-means se
  s2_trt_2018 <- var(trt_2018)
  s2_trt_2017 <- var(trt_2017)
  s2_ctrl_2018 <- var(ctrl_2018)   # pooled control variance at post
  s2_ctrl_2017 <- var(ctrl_2017)   # pooled control variance at pre
  
  # pooled se
  se_diff_in_means <- sqrt(
    s2_trt_2018/length(trt_2018)+
      s2_trt_2017/length(trt_2017)+
      s2_ctrl_2018/length(ctrl_2018)+
      s2_ctrl_2017/length(ctrl_2017)
    )
  
  pl_hldr <- rep(99, 38)
  cov_diff_in_means <- c(pl_hldr, true_tau >= diff_in_means - 1.96*se_diff_in_means & true_tau <= diff_in_means + 1.96*se_diff_in_means)
  diff_in_means <- c(pl_hldr, diff_in_means)
  se_diff_in_means <- c(pl_hldr, se_diff_in_means)
  names(se_diff_in_means) <- c(1979:2016, 2018)
  names(diff_in_means) <- c(1979:2016, 2018)
  names(cov_diff_in_means) <- c(1979:2016, 2018)
  
  # run two-group, multiple pre-period event study
  time_start <- Sys.time()
  lm_ES <- lm(Y_i ~ D_i+T_f+D_i:T_f, data=rcs_data)
  
  # get heteroskedasticity-robust (Eicker–Huber–White) SE, but assumes independent sampling across i so no clustering
  SE <- vcovHC(lm_ES, type="HC0")
  time_end <- Sys.time()
  time_elapsed_pt <- rep(as.numeric(difftime(time_end, time_start, units = "secs")), 39)
  names(time_elapsed_pt) <- c(1979:2016, 2018)
  
  # organize results
  beta_idx <- grepl("^D_i:T_f", names(coef(lm_ES)))
  tau_pt <- coef(lm_ES)[beta_idx]
  se_pt <- sqrt(diag(SE))[beta_idx]
  coverage_pt <- as.numeric(true_tau >= tau_pt - 1.96*se_pt & true_tau <= tau_pt + 1.96*se_pt)
  
  names(tau_pt) <- c(1979:2016, 2018)
  names(se_pt) <- c(1979:2016, 2018)
  names(coverage_pt) <- c(1979:2016, 2018)
  
  # collect results
  result_pt <- bind_rows(
    result_pt,
    c(as.list(tau_pt), list(type = "est", DGP = 1, niter=iter)),
    c(as.list(se_pt), list(type = "se", DGP = 1, niter=iter)),
    c(as.list(coverage_pt), list(type = "coverage", DGP = 1, niter=iter)),
    c(as.list(time_elapsed_pt), list(type = "time", DGP = 1, niter=iter)),
    c(as.list(diff_in_means), list(type = "est_dim", DGP = 1, niter=iter)),
    c(as.list(se_diff_in_means), list(type = "se_dim", DGP = 1, niter=iter)),
    c(as.list(cov_diff_in_means), list(type = "coverage_dim", DGP = 1, niter=iter))
  )
  ## END: PT-DID ------------------------------------
  
  
  ## START: SDID ------------------------------------
  # compute sample means
  Y_means <- rcs_data %>%
    group_by(G_i, T_i) %>%
    summarise(wage=mean(Y_i), .groups = "drop")
  # treatment status
  Y_means$assignment <- Y_means$G_i>50 & Y_means$T_i==2018
  # prep data into matrix format
  SDID_setup <- panel.matrices(as.data.frame(Y_means),
                               unit="G_i",
                               time="T_i",
                               treatment='assignment',
                               outcome='wage')
  time_start <- Sys.time()
  
  method_sdid <- synthdid_estimate(SDID_setup$Y, 
                                   SDID_setup$N0,
                                   SDID_setup$T0)
  tau_sdid <- summary(method_sdid)$estimate
  se_sdid <- sqrt(vcov(method_sdid, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sdid <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sdid <- as.numeric(true_tau >= tau_sdid - 1.96*se_sdid & true_tau <= tau_sdid + 1.96*se_sdid)
  
  # collect results
  result_sdid <- bind_rows(
    result_sdid,
    list("2018" = tau_sdid, type = "est", DGP = 1, niter=iter),
    list("2018" = se_sdid, type = "se", DGP = 1, niter=iter),
    list("2018" = coverage_sdid, type = "coverage", DGP = 1, niter=iter),
    list("2018" = time_elapsed_sdid, type = "time", DGP = 1, niter=iter)
  )
  ## END: SDID --------------------------------------
  
  
  ## START: SC ------------------------------------
  time_start <- Sys.time()
  
  method_sc <- sc_estimate(SDID_setup$Y,
                           SDID_setup$N0,
                           SDID_setup$T0)
  tau_sc <- summary(method_sc)$estimate
  se_sc <- sqrt(vcov(method_sc, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sc <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sc <- as.numeric(true_tau >= tau_sc - 1.96*se_sc & true_tau <= tau_sc + 1.96*se_sc)
  
  # collect results
  result_sc <- bind_rows(
    result_sc,
    list("2018" = tau_sc, type = "est", DGP = 1, niter=iter),
    list("2018" = se_sc, type = "se", DGP = 1, niter=iter),
    list("2018" = coverage_sc, type = "coverage", DGP = 1, niter=iter),
    list("2018" = time_elapsed_sc, type = "time", DGP = 1, niter=iter)
  )
  ## END: SC ------------------------------------
  
  
  ## START: SPT ------------------------------------
  time_start <- Sys.time()
  n <- length(rcs_data$Y_i)
  
  method_spt <- tau_set(outcome=rcs_data$Y_i,
                        unit_id=rcs_data$G_i,
                        time_id=rcs_data$T_i,
                        trt_unit_id=which(trt_idx==1)+50,
                        cand_tau_size=cand_tau_size,
                        num_bstp_rep=num_bstp_rep,
                        cand_tau=cand_tau_PT)
  
  time_end <- Sys.time()
  time_elapsed_spt <- as.numeric(difftime(time_end, time_start, units = "secs"))
  cs_tau <- method_spt$cand_tau[method_spt$rej==0]
  
  # for future iterations, candidate values of tau for the next iteration are based on enlarging current iteration's cs_tau
  # double the distance from the midpoint (symmetric expansion)
  a <- min(cs_tau)
  b <- max(cs_tau)
  mid <- (a+b)/2
  len <- b-a
  a <- mid-len
  b <- mid+len
  cand_tau_PT <- a+(b-a)*
    qbeta(seq(0, 1, length.out = cand_tau_size), shape1 = 2, shape2 = 2)
  
  # collect results
  result_spt <- bind_rows(
    result_spt,
    list("2018" = min(cs_tau), type = "cs_lb", DGP = 1, niter=iter),
    list("2018" = max(cs_tau), type = "cs_ub", DGP = 1, niter=iter),
    list("2018" = as.numeric(true_tau <= max(cs_tau) & true_tau >= min(cs_tau)), type = "coverage", DGP = 1, niter=iter),
    list("2018" = time_elapsed_spt, type = "time", DGP = 1, niter=iter)
  )
  frac_rej_PT <- append(frac_rej_PT, 1-length(cs_tau)/cand_tau_size)
  ## END: SPT ------------------------------------
  ### END: DGP-1 treated unit is one state, PT is satisfied ----
  
  
  ## START: DGP-2 treated unit is one state, PT is satisfied for pre-trends but not for post-trend ----
  # compute population means
  pop_mean_mat <- cell_stats %>% mutate(W=state_num %in% which(trt_idx==1) & year==2018)
  pop_mean_mat <- panel.matrices(data.frame(pop_mean_mat), outcome="wage", treatment="W")$Y
  pop_trend_mat <- pop_mean_mat[, 2:40] - pop_mean_mat[, 1:39]
  
  A_pre <- pop_trend_mat[1:49, 1:38] # 49 by 38
  A_pre <- A_pre[order(as.numeric(row.names(A_pre))), ] # sort control indices in increasing order, as in cell_means
  A_pre <- t(A_pre) # 38 x 49
  b_pre <- matrix(trt_unit$trt_trends[1:38]) # the simulated pre-trends for the treated unit, 1 by 38
  a_post <- pop_trend_mat[1:49, 39, drop=FALSE] # 49 x 1
  a_post <- a_post[order(as.numeric(row.names(a_post))), , drop=FALSE]
  PT1_ctfact_trend <- trt_unit$trt_trends[39] # counterfactual trend under DGP-PT1
  
  ## SANITY CHECK: Apre pi_0 = b_pre
  # stopifnot(max(abs(A_pre %*% pi_0 - b_pre)) < 1e-10,
  #           abs(sum(pi_0) - 1) < 1e-10,
  #           all(pi_0 >= -1e-10))
  
  K0 <- ncol(A_pre) # number of control units
  T0 <- nrow(A_pre)+1 # number of pre-periods
  
  for (trial in 1:20000){
    # try a different a_post that has more dispersion so we get a meaningful violation of post-trend
    a_post_alt <- a_post + rnorm(length(a_post), 0, 4*se_pt[length(se_pt)]) 
    # set violation of post-trend to be 2x se from event study
    get_w_alt <- get_alt_omega(A_pre, b_pre, a_post_alt, pi_0,
                               delta_target=2*se_pt[length(se_pt)])
    if(!is.null(get_w_alt$w_alt)) break
  }
  
  a_post_alt <- get_w_alt$a_post_alt
  names(a_post_alt) <- rownames(a_post)
  w_alt <- get_w_alt$w_alt
  names(w_alt) <- rownames(a_post)
  violation_noPT <- append(violation_noPT, get_w_alt$violation)
  
  # ## SANITY CHECK: A_pre w_alt = b_pre
  # stopifnot(max(abs(A_pre %*% w_alt - b_pre)) < 1e-6,
  #           abs(sum(w_alt) - 1) < 1e-6,
  #           all(w_alt >= -1e-6))
  
  # generate new post-treatment outcome for controls
  ctrl_mean_T0 <- pop_mean_mat[1:49, 39, drop=FALSE]
  ctrl_mean_T0 <- ctrl_mean_T0[order(as.numeric(row.names(ctrl_mean_T0))), , drop=FALSE] # 49 by 1, T0=2017 outcomes
  ctrl_mean_post_noPT <- ctrl_mean_T0 + a_post_alt
  
  # generate repeated cross-sections for the control units with the updated post means
  rcs_data_ctrl_noPT <- map2_dfr(which(trt_idx==0), n_k0, \(k, nk)
                                 map_dfr(1979:2018, \(t)
                                         gen_micro_data_control(k, t, nk, PT=FALSE)))
  
  # generate repeated cross-sections for the treated unit, now using the w_alt weighting so that PT fails in the post-period
  trt_unit_noPT <- gen_trt_unit(trt_idx, pi_0=w_alt, PT=FALSE)
  rcs_data_trt_noPT <- map2_dfr(which(trt_idx==1), n_k1, \(k, nk)
                                map_dfr(1979:2018, \(t)
                                        gen_micro_data_trt(trt_info=trt_unit_noPT,
                                                           k=k, t=t, n=nk,
                                                           true_tau=true_tau)))
  
  # combine into one data frame
  rcs_data_noPT <- bind_rows(rcs_data_ctrl_noPT, rcs_data_trt_noPT)
  
  # set T_0=2017 to be the reference period
  rcs_data_noPT$T_f <- factor(rcs_data_noPT$T_i, 
                              levels = sort(unique(rcs_data_noPT$T_i)))
  rcs_data_noPT$T_f <- relevel(rcs_data_noPT$T_f, ref = "2017")  
  
  ## START: PT-DID ------------------------------------
  # double check 2018 coeff matches difference-in-means estimator
  trt_2018 <- subset(rcs_data_noPT, T_i==2018 & D_i==1)$Y_i
  trt_2017 <- subset(rcs_data_noPT, T_i==2017 & D_i==1)$Y_i
  ctrl_2018 <- subset(rcs_data_noPT, T_i==2018 & D_i==0)$Y_i
  ctrl_2017 <- subset(rcs_data_noPT, T_i==2017 & D_i==0)$Y_i
  diff_in_means <- mean(trt_2018)-mean(trt_2017)-(mean(ctrl_2018)-mean(ctrl_2017))
  
  # difference-in-means se
  s2_trt_2018 <- var(trt_2018)
  s2_trt_2017 <- var(trt_2017)
  s2_ctrl_2018 <- var(ctrl_2018)   # pooled control variance at post
  s2_ctrl_2017 <- var(ctrl_2017)   # pooled control variance at pre
  
  # pooled se
  se_diff_in_means <- sqrt(s2_trt_2018/length(trt_2018)+
                             s2_trt_2017/length(trt_2017)+ 
                             s2_ctrl_2018/length(ctrl_2018)+ 
                             s2_ctrl_2017/length(ctrl_2017))
  
  cov_diff_in_means <- c(pl_hldr, true_tau >= diff_in_means - 1.96*se_diff_in_means & true_tau <= diff_in_means + 1.96*se_diff_in_means)
  diff_in_means <- c(pl_hldr, diff_in_means)
  se_diff_in_means <- c(pl_hldr, se_diff_in_means)
  names(diff_in_means) <- c(1979:2016, 2018)
  names(se_diff_in_means) <- c(1979:2016, 2018)
  names(cov_diff_in_means) <- c(1979:2016, 2018)
  
  # run two-group, multiple pre-period event study
  time_start <- Sys.time()
  lm_ES <- lm(Y_i ~ D_i+T_f+D_i:T_f, data=rcs_data_noPT)
  
  # get heteroskedasticity-robust (Eicker–Huber–White) SE, but assumes independent sampling across i so no clustering
  SE <- vcovHC(lm_ES, type="HC0")
  time_end <- Sys.time()
  time_elapsed_pt <- rep(as.numeric(difftime(time_end, time_start, units = "secs")), 39)
  names(time_elapsed_pt) <- c(1979:2016, 2018)
  
  # organize results
  beta_idx <- grepl("^D_i:T_f", names(coef(lm_ES)))
  tau_pt <- coef(lm_ES)[beta_idx]
  se_pt <- sqrt(diag(SE))[beta_idx]
  coverage_pt <- as.numeric(true_tau >= tau_pt - 1.96*se_pt & true_tau <= tau_pt + 1.96*se_pt)
  
  names(tau_pt) <- c(1979:2016, 2018)
  names(se_pt) <- c(1979:2016, 2018)
  names(coverage_pt) <- c(1979:2016, 2018)
  
  # collect results
  result_pt <- bind_rows(
    result_pt,
    c(as.list(tau_pt), list(type = "est", DGP = 2, niter=iter)),
    c(as.list(se_pt), list(type = "se", DGP = 2, niter=iter)),
    c(as.list(coverage_pt), list(type = "coverage", DGP = 2, niter=iter)),
    c(as.list(time_elapsed_pt), list(type = "time", DGP = 2, niter=iter)),
    c(as.list(diff_in_means), list(type = "est_dim", DGP = 2, niter=iter)),
    c(as.list(se_diff_in_means), list(type = "se_dim", DGP = 2, niter=iter)),
    c(as.list(cov_diff_in_means), list(type = "coverage_dim", DGP = 2, niter=iter))
  )
  ## END: PT-DID ------------------------------------
  
  
  ## START: SDID ------------------------------------
  # compute sample means
  Y_means <- rcs_data_noPT %>%
    group_by(G_i, T_i) %>%
    summarise(wage=mean(Y_i), .groups = "drop")
  # treatment status
  Y_means$assignment <- Y_means$G_i>50 & Y_means$T_i==2018
  # prep data into matrix format
  SDID_setup <- panel.matrices(as.data.frame(Y_means),
                               unit="G_i",
                               time="T_i",
                               treatment='assignment',
                               outcome='wage')
  time_start <- Sys.time()
  
  method_sdid <- synthdid_estimate(SDID_setup$Y, 
                                   SDID_setup$N0,
                                   SDID_setup$T0)
  tau_sdid <- summary(method_sdid)$estimate
  se_sdid <- sqrt(vcov(method_sdid, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sdid <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sdid <- as.numeric(true_tau >= tau_sdid - 1.96*se_sdid & true_tau <= tau_sdid + 1.96*se_sdid)
  
  # collect results
  result_sdid <- bind_rows(
    result_sdid,
    list("2018" = tau_sdid, type = "est", DGP = 2, niter=iter),
    list("2018" = se_sdid, type = "se", DGP = 2, niter=iter),
    list("2018" = coverage_sdid, type = "coverage", DGP = 2, niter=iter),
    list("2018" = time_elapsed_sdid, type = "time", DGP = 2, niter=iter)
  )
  ## END: SDID --------------------------------------
  
  
  ## START: SC ------------------------------------
  time_start <- Sys.time()
  
  method_sc <- sc_estimate(SDID_setup$Y,
                           SDID_setup$N0,
                           SDID_setup$T0)
  tau_sc <- summary(method_sc)$estimate
  se_sc <- sqrt(vcov(method_sc, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sc <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sc <- as.numeric(true_tau >= tau_sc - 1.96*se_sc & true_tau <= tau_sc + 1.96*se_sc)
  
  # collect results
  result_sc <- bind_rows(
    result_sc,
    list("2018" = tau_sc, type = "est", DGP = 2, niter=iter),
    list("2018" = se_sc, type = "se", DGP = 2, niter=iter),
    list("2018" = coverage_sc, type = "coverage", DGP = 2, niter=iter),
    list("2018" = time_elapsed_sc, type = "time", DGP = 2, niter=iter)
  )
  ## END: SC ------------------------------------
  
  
  ## START: SPT ------------------------------------
  time_start <- Sys.time()
  n <- length(rcs_data_noPT$Y_i)
  method_spt <- tau_set(outcome=rcs_data_noPT$Y_i,
                        unit_id=rcs_data_noPT$G_i,
                        time_id=rcs_data_noPT$T_i,
                        trt_unit_id=which(trt_idx==1)+50,
                        cand_tau_size=cand_tau_size,
                        num_bstp_rep=num_bstp_rep,
                        cand_tau=cand_tau_noPT)
  
  time_end <- Sys.time()
  time_elapsed_spt <- as.numeric(difftime(time_end, time_start, units = "secs"))
  cs_tau <- method_spt$cand_tau[method_spt$rej==0]
  
  # for future iterations, candidate values of tau for the next iteration are based on enlarging current iteration's cs_tau
  # double the distance from the midpoint (symmetric expansion)
  a <- min(cs_tau)
  b <- max(cs_tau)
  mid <- (a+b)/2
  len <- b-a
  a <- mid-len
  b <- mid+len
  cand_tau_noPT <- a+(b-a)*
    qbeta(seq(0, 1, length.out = cand_tau_size), shape1 = 2, shape2 = 2)
  
  # collect results
  result_spt <- bind_rows(
    result_spt,
    list("2018" = min(cs_tau), type = "cs_lb", DGP = 2, niter=iter),
    list("2018" = max(cs_tau), type = "cs_ub", DGP = 2, niter=iter),
    list("2018" = as.numeric(true_tau <= max(cs_tau) & true_tau >= min(cs_tau)), type = "coverage", DGP = 2, niter=iter),
    list("2018" = time_elapsed_spt, type = "time", DGP = 2, niter=iter)
  )
  frac_rej_noPT <- append(frac_rej_noPT, 1-length(cs_tau)/cand_tau_size)
  ## END: SPT ------------------------------------
  ### END: DGP-2 treated unit is one state, PT is satisfied for pre-trends but not for post-trend ----
  
  end <- Sys.time()
  
  print(paste0("DGP-1: PT holds --------------"))
  print(paste0("PT-DID avg bias: ", mean(abs(true_tau-subset(result_pt, type=="est" & DGP==1)$`2018`))))
  print(paste0("PT-DID avg bias (DiM): ", mean(abs(true_tau-subset(result_pt, type=="est_dim" & DGP==1)$`2018`))))
  print(paste0("PT-DID avg CI length: ", 1.96*2*mean(subset(result_pt, type=="se" & DGP==1)$`2018`)))
  print(paste0("PT-DID avg CI length (DiM): ", 1.96*2*mean(subset(result_pt, type=="se_dim" & DGP==1)$`2018`)))
  print(paste0("PT-DID coverage rate: ", mean(subset(result_pt, type=="coverage" & DGP==1)$`2018`)))
  print(paste0("PT-DID coverage rate (DiM): ", mean(subset(result_pt, type=="coverage_dim" & DGP==1)$`2018`)))
  print(paste0("PT-DID avg time: ", mean(subset(result_pt, type=="time" & DGP==1)$`2018`)))
  print(paste0("-"))
  print(paste0("SDID avg bias: ", mean(abs(true_tau-subset(result_sdid, type=="est" & DGP==1)$`2018`))))
  print(paste0("SDID avg CI length: ", 1.96*2*mean(subset(result_sdid, type=="se" & DGP==1)$`2018`)))
  print(paste0("SDID coverage rate: ", mean(subset(result_sdid, type=="coverage" & DGP==1)$`2018`)))
  print(paste0("SDID avg time: ", mean(subset(result_sdid, type=="time" & DGP==1)$`2018`)))
  print(paste0("-"))
  print(paste0("SC avg bias: ", mean(abs(true_tau-subset(result_sc, type=="est" & DGP==1)$`2018`))))
  print(paste0("SC avg CI length: ", 1.96*2*mean(subset(result_sc, type=="se" & DGP==1)$`2018`)))
  print(paste0("SC coverage rate: ", mean(subset(result_sc, type=="coverage" & DGP==1)$`2018`)))
  print(paste0("SC avg time: ", mean(subset(result_sc, type=="time" & DGP==1)$`2018`)))
  print(paste0("-"))
  print(paste0("SPT avg CS length: ", mean(subset(result_spt, type=="cs_ub" & DGP==1)$`2018`-subset(result_spt, type=="cs_lb" & DGP==1)$`2018`)))
  print(paste0("SPT coverage rate: ", mean(subset(result_spt, type=="coverage" & DGP==1)$`2018`)))
  print(paste0("SPT avg time: ", mean(subset(result_spt, type=="time" & DGP==1)$`2018`)))
  print(paste0("frac cand_tau rej: ", mean(frac_rej_PT)))
  print(paste0("DGP-2: Post PT fails --------------"))
  print(paste0("avg dev from parallel posttrend: ", mean(violation_noPT)))
  print(paste0("PT-DID avg bias: ", mean(abs(true_tau-subset(result_pt, type=="est" & DGP==2)$`2018`))))
  print(paste0("PT-DID avg bias (DiM): ", mean(abs(true_tau-subset(result_pt, type=="est_dim" & DGP==2)$`2018`))))
  print(paste0("PT-DID avg CI length: ", 1.96*2*mean(subset(result_pt, type=="se" & DGP==2)$`2018`)))
  print(paste0("PT-DID avg CI length (DiM): ", 1.96*2*mean(subset(result_pt, type=="se_dim" & DGP==2)$`2018`)))
  print(paste0("PT-DID coverage rate: ", mean(subset(result_pt, type=="coverage" & DGP==2)$`2018`)))
  print(paste0("PT-DID coverage rate (DiM): ", mean(subset(result_pt, type=="coverage_dim" & DGP==2)$`2018`)))
  print(paste0("PT-DID avg time: ", mean(subset(result_pt, type=="time" & DGP==2)$`2018`)))
  print(paste0("-"))
  print(paste0("SDID avg bias: ", mean(abs(true_tau-subset(result_sdid, type=="est" & DGP==2)$`2018`))))
  print(paste0("SDID avg CI length: ", 1.96*2*mean(subset(result_sdid, type=="se" & DGP==2)$`2018`)))
  print(paste0("SDID coverage rate: ", mean(subset(result_sdid, type=="coverage" & DGP==2)$`2018`)))
  print(paste0("SDID avg time: ", mean(subset(result_sdid, type=="time" & DGP==2)$`2018`)))
  print(paste0("-"))
  print(paste0("SC avg bias: ", mean(abs(true_tau-subset(result_sc, type=="est" & DGP==2)$`2018`))))
  print(paste0("SC avg CI length: ", 1.96*2*mean(subset(result_sc, type=="se" & DGP==2)$`2018`)))
  print(paste0("SC coverage rate: ", mean(subset(result_sc, type=="coverage" & DGP==2)$`2018`)))
  print(paste0("SC avg time: ", mean(subset(result_sc, type=="time" & DGP==2)$`2018`)))
  print(paste0("-"))
  print(paste0("SPT avg CS length: ", mean(subset(result_spt, type=="cs_ub" & DGP==2)$`2018`-subset(result_spt, type=="cs_lb" & DGP==2)$`2018`)))
  print(paste0("SPT coverage rate: ", mean(subset(result_spt, type=="coverage" & DGP==2)$`2018`)))
  print(paste0("SPT avg time: ", mean(subset(result_spt, type=="time" & DGP==2)$`2018`)))
  print(paste0("frac cand_tau rej: ", mean(frac_rej_noPT)))
  
  print(end-start)
  
  write.csv(result_spt, "../result/result_spt_DGP12.csv", row.names=FALSE)
  write.csv(result_pt, "../result/result_pt_DGP12.csv", row.names=FALSE)
  write.csv(result_sc, "../result/result_sc_DGP12.csv", row.names=FALSE)
  write.csv(result_sdid, "../result/result_sdid_DGP12.csv", row.names=FALSE)
}
