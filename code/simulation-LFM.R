library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(sandwich)
library(lmtest)
library(broom)
library(osqp)
library(CVXR)
library(Matrix)
library(mvtnorm)
library(pracma)
library(synthdid)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
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
  #### ------ START OF MODIFICATION ------
  E <- mvtnorm::rmvnorm(N, sigma = Sigma, method = "chol")
  Y =  L + E 
  rownames(Y) <- 1:50
  # return(list(Y=Y[order(assignment), ], N0=N0, T0=T0))
  return(list(Y=Y[order(assignment), ], N0=N0, T0=T0, assignment=assignment))
  #### ------ END OF MODIFICATION ------ 
  # order units by treatment, so treated units are at the bottom of our data matrix in accordance with our convention
}

# fit a rank-4 factor matrix L to the outcome; decompose L = F + M; estimate error AR2 covariance E; estimate assignment probability
factor_model=estimate_dgp(Y.wage, w.minwage, rank=4)
# get the covariance matrix for E
Sigma <- factor_model$Sigma
# get the factor component L
L <- factor_model$`F`+factor_model$M
rownames(L) <- 1:50
colnames(L) <- 1979:2018

# NO rows of L can be expressed as a convex combination of the other 49 rows under the factor component L. So we will need to impute one row in the DGP-SDID to be the convex combination of the remaining rows so that SPT doesn't return an empty CS (which it will, if we don't do this).
K <- nrow(L)
convex_state <- c()
for (k in 1:K){
  target <- L[k, ]                   
  others <- L[-k, , drop = FALSE] 
  
  w <- Variable(K-1)  # weights on other rows
  
  constraints <- list(
    t(others) %*% w == target,  # convex combination equality
    sum(w) == 1,  w >= 0                
  )
  
  prob <- Problem(Minimize(0), constraints)
  res <- suppressWarnings(
    solve(prob, solver="ECOS",
          ecos.opts=list(reltol=1e-6, 
                         abstol=1e-6, 
                         feastol=1e-6))
    )
  if (res$status == "optimal"){
    convex_state <- append(convex_state, 1)
  } else{
    convex_state <- append(convex_state, 0)
  }
}

print(paste0("There are ", sum(convex_state), " states that can be expressed as a convex combination of the other 49 states under the factor component L."))

# simulation configuration
true_tau <- 0 # null effect
num_MC <- 1000
cand_tau_size <- 500
num_bstp_rep <- 1000
cand_tau_SDID <- NULL
cand_tau_noSDID <- NULL
frac_rej_SDID <- c()
frac_rej_noSDID <- c()
violation_noPT <- c()

# collect results; type in c("est", "se", "coverage", "cs_lb", "cs_ub", "time", "est_lb", "est_ub"); DGP in c(3, 4)
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
  set.seed(iter)
  print(paste0("iter:  ", iter))
  start <- Sys.time()
  ### START: DGP-3 treated unit is one state, many weights are valid pre and post ----
  # simulate Y = L + true_tau*W + E for SDID
  factor_model_1tr <- simulate_dgp(factor_model, N1=1, T1=1)
  
  # get the treated unit index
  trt_idx <- factor_model_1tr$assignment
  
  # randomly draw a convex weight
  convex_w <- rexp(49)
  convex_w <- convex_w/sum(convex_w)
  
  # generate a fake treated factor equal to the convex comb of the control factors
  L <- factor_model$`F`+factor_model$M
  L[which(trt_idx==1), ] <- t(L[which(trt_idx==0), , drop = FALSE])%*%matrix(convex_w)
  
  trt_unit <- gen_trt_unit(trt_idx=trt_idx)
  pi_0 <- trt_unit$pi_0
  pi_1 <- trt_unit$pi_1
  
  # generate micro-level
  n_1 <- sum(subset(cell_stats, year==2018 & state_num %in% which(trt_idx==1))$n_kt) # total number of treated individuals
  n_0 <- sum(subset(cell_stats, year==2018 & state_num %in% which(trt_idx==0))$n_kt) # total number of control individuals
  n_k0 <- ceiling(n_0*pi_0) # control individuals sampled proportional to pi_0
  n_k1 <- ceiling(n_1*pi_1) # treated individuals sampled proportional to pi_1
  
  # panel for the control units
  pnl_data_ctrl <- map2_dfr(which(trt_idx == 0), n_k0, \(k, nk)
                            gen_micro_data_control(k=k, n=nk, PT=FALSE, 
                                                   SDID=TRUE, Sigma=Sigma))
  
  # panel for the treated unit
  pnl_data_trt <- map2_dfr(which(trt_idx==1), n_k1, \(k, nk)
                           gen_micro_data_trt(k=k, n=nk,
                                              true_tau=true_tau,
                                              SDID=TRUE, Sigma=Sigma))
  
  # combine into one data frame
  pnl_data <- bind_rows(pnl_data_ctrl, pnl_data_trt)
  
  # set T_0=2017 to be the reference period
  pnl_data$T_f <- factor(pnl_data$T_i, 
                         levels = sort(unique(pnl_data$T_i)))
  pnl_data$T_f <- relevel(pnl_data$T_f, ref = "2017")  
  
  ## START: PT-DID ------------------------------------
  # get first order differences Y_{i,2018}-Y_{i,2017} within state-year cell
  FoD <- pnl_data %>%
    filter(T_i %in% c(2017, 2018)) %>%
    arrange(G_i, person_id, T_i) %>%
    group_by(G_i, person_id, D_i) %>%
    filter(n() == 2) %>%
    summarise(dy = diff(Y_i), .groups = "drop") %>% 
    group_by(D_i) %>%
    summarise(
      n = n(), diff_mean = mean(dy), s2 = var(dy),
      .groups = "drop"
    )
  
  # difference-in-means estimator & SE
  diff_in_means <- subset(FoD, D_i==1)$diff_mean-subset(FoD, D_i==0)$diff_mean
  se_diff_in_means <- sqrt(subset(FoD, D_i==1)$s2*(subset(FoD, D_i==1)$n-1)/(subset(FoD, D_i==1)$n)^2+subset(FoD, D_i==0)$s2*(subset(FoD, D_i==0)$n-1)/(subset(FoD, D_i==0)$n)^2) # var() uses denominator (n_g-1); i'd like n_g for comparison with HC0
  
  pl_hldr <- rep(99, 38)
  cov_diff_in_means <- c(pl_hldr, true_tau >= diff_in_means - 1.96*se_diff_in_means & true_tau <= diff_in_means + 1.96*se_diff_in_means)
  diff_in_means <- c(pl_hldr, diff_in_means)
  se_diff_in_means <- c(pl_hldr, se_diff_in_means)
  names(se_diff_in_means) <- c(1979:2016, 2018)
  names(diff_in_means) <- c(1979:2016, 2018)
  names(cov_diff_in_means) <- c(1979:2016, 2018)
  
  # run two-group, multiple pre-period event study
  time_start <- Sys.time()
  lm_ES <- lm(Y_i ~ D_i+T_f+D_i:T_f, data=pnl_data)
  
  # get heteroskedasticity-robust (Eicker–Huber–White) SE clustered at the person level for panel data (independent across i and states, with serial correlation across time for each i)
  SE <- vcovCL(lm_ES, cluster = ~person_id, type = "HC0")
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
    c(as.list(tau_pt), list(type = "est", DGP = 3, niter=iter)),
    c(as.list(se_pt), list(type = "se", DGP = 3, niter=iter)),
    c(as.list(coverage_pt), list(type = "coverage", DGP = 3, niter=iter)),
    c(as.list(time_elapsed_pt), list(type = "time", DGP = 3, niter=iter)),
    c(as.list(diff_in_means), list(type = "est_dim", DGP = 3, niter=iter)),
    c(as.list(se_diff_in_means), list(type = "se_dim", DGP = 3, niter=iter)),
    c(as.list(cov_diff_in_means), list(type = "coverage_dim", DGP = 3, niter=iter))
  )
  ## END: PT-DID ------------------------------------
  
  ## START: SDID ------------------------------------
  time_start <- Sys.time()
  
  method_sdid <- synthdid_estimate(factor_model_1tr$Y, 
                                   factor_model_1tr$N0,
                                   factor_model_1tr$T0)
  tau_sdid <- summary(method_sdid)$estimate
  se_sdid <- sqrt(vcov(method_sdid, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sdid <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sdid <- as.numeric(true_tau >= tau_sdid - 1.96*se_sdid & true_tau <= tau_sdid + 1.96*se_sdid)
  
  # collect results
  result_sdid <- bind_rows(
    result_sdid,
    list("2018" = tau_sdid, type = "est", DGP = 3, niter=iter),
    list("2018" = se_sdid, type = "se", DGP = 3, niter=iter),
    list("2018" = coverage_sdid, type = "coverage", DGP = 3, niter=iter),
    list("2018" = time_elapsed_sdid, type = "time", DGP = 3, niter=iter)
  )
  ## END: SDID --------------------------------------
  
  ## START: SC ------------------------------------
  time_start <- Sys.time()
  
  method_sc <- sc_estimate(factor_model_1tr$Y, 
                           factor_model_1tr$N0,
                           factor_model_1tr$T0)
  tau_sc <- summary(method_sc)$estimate
  se_sc <- sqrt(vcov(method_sc, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sc <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sc <- as.numeric(true_tau >= tau_sc - 1.96*se_sc & true_tau <= tau_sc + 1.96*se_sc)
  
  # collect results
  result_sc <- bind_rows(
    result_sc,
    list("2018" = tau_sc, type = "est", DGP = 3, niter=iter),
    list("2018" = se_sc, type = "se", DGP = 3, niter=iter),
    list("2018" = coverage_sc, type = "coverage", DGP = 3, niter=iter),
    list("2018" = time_elapsed_sc, type = "time", DGP = 3, niter=iter)
  )
  ## END: SC --------------------------------------
  
  
  ## START: SPT ------------------------------------
  time_start <- Sys.time()
  n <- length(pnl_data$Y_i)
  
  method_spt <- tau_set(outcome=pnl_data$Y_i,
                        unit_id=pnl_data$G_i,
                        time_id=pnl_data$T_i,
                        trt_unit_id=which(trt_idx==1)+50,
                        person_id=pnl_data$person_id,
                        cand_tau_size=cand_tau_size,
                        num_bstp_rep=num_bstp_rep,
                        cand_tau=cand_tau_SDID)
  
  time_end <- Sys.time()
  time_elapsed_spt <- as.numeric(difftime(time_end, time_start, units = "secs"))
  cs_tau <- method_spt$cand_tau[method_spt$rej==0]
  
  # for future iterations, candidate values of tau for the next iteration are based on enlarging current iteration's cs_tau, if nonempty
  if (length(cs_tau)>0){
    # double the distance from the midpoint (symmetric expansion)
    a <- min(cs_tau)
    b <- max(cs_tau)
    mid <- (a+b)/2
    len <- b-a
    a <- mid-len
    b <- mid+len
    cand_tau_SDID <- a+(b-a)*
      qbeta(seq(0, 1, length.out = cand_tau_size), shape1 = 2, shape2 = 2)
  }
  
  # collect results
  result_spt <- bind_rows(
    result_spt,
    list("2018" = ifelse(length(cs_tau)>0, min(cs_tau), -99), type = "cs_lb", DGP = 3, niter=iter),
    list("2018" = ifelse(length(cs_tau)>0, max(cs_tau), -99), type = "cs_ub", DGP = 3, niter=iter),
    list("2018" = ifelse(length(cs_tau)>0, as.numeric(true_tau <= max(cs_tau) & true_tau >= min(cs_tau)), 0), type = "coverage", DGP = 3, niter=iter),
    list("2018" = time_elapsed_spt, type = "time", DGP = 3, niter=iter)
  )
  frac_rej_SDID <- append(frac_rej_SDID, 1-length(cs_tau)/cand_tau_size)
  ## END: SPT ------------------------------------
  ### END: DGP-3 treated unit is one state, many weights are valid pre and post ----
  
  
  
  ## START: DGP-4 treated unit is one state, many weights are valid pre but only a sparse weight is valid post----
  trt_i <- which(trt_idx==1) # trt unit's index
  controls <- setdiff(1:50, trt_i) # control indices
  D <- 3 # population true donors count
  donors <- sample(controls, D) # randomly sample 3 donors
  clones <- setdiff(controls, donors)  # the rest control units (clones) are irrelevant to the treated unit
  C <- length(clones)
  
  L <- factor_model$`F`+factor_model$M
  
  # treated equals donors’ average in ALL periods
  L[trt_i, ] <- colMeans(L[donors, , drop=FALSE])
  
  # counterfactual given by the donors
  ctfactl <- L[trt_i, 40]
  
  # reset all clones' pre-trt factors to the treated unit's factors: all convex weights are valid
  L[clones, 1:39] <- matrix(
    rep(L[trt_i, 1:39], times = length(clones)),
    nrow = length(clones), byrow = TRUE
  )
  
  # in post period, the other irrelevant controls/clones are systematically different from the counterfactual by -5
  L[clones, 40] <- ctfactl - 5
  
  # regenerate panel for control units
  pnl_data_ctrl_noSDID <- map2_dfr(which(trt_idx == 0), n_k0, \(k, nk)
                            gen_micro_data_control(k=k, n=nk, PT=FALSE, 
                                                   SDID=TRUE,
                                                   Sigma=Sigma))
  
  # regenerate panel for the treated unit
  pnl_data_trt_noSDID <- map2_dfr(which(trt_idx==1), n_k1, \(k, nk)
                                  gen_micro_data_trt(k=k, n=nk,
                                                     true_tau=true_tau,
                                                     SDID=TRUE, 
                                                     Sigma=Sigma))
  
  # combine into one data frame
  pnl_data_noSDID <- bind_rows(pnl_data_ctrl_noSDID, pnl_data_trt_noSDID)
  
  # set T_0=2017 to be the reference period
  pnl_data_noSDID$T_f <- factor(pnl_data_noSDID$T_i, 
                         levels = sort(unique(pnl_data_noSDID$T_i)))
  pnl_data_noSDID$T_f <- relevel(pnl_data_noSDID$T_f, ref = "2017")  
  
  ## START: PT-DID ------------------------------------
  # get first order differences Y_{i,2018}-Y_{i,2017} within state-year cell
  FoD <- pnl_data_noSDID %>%
    filter(T_i %in% c(2017, 2018)) %>%
    arrange(G_i, person_id, T_i) %>%
    group_by(G_i, person_id, D_i) %>%
    filter(n() == 2) %>%
    summarise(dy = diff(Y_i), .groups = "drop") %>% 
    group_by(D_i) %>%
    summarise(
      n = n(), diff_mean = mean(dy), s2 = var(dy),
      .groups = "drop"
    )
  
  # difference-in-means estimator & SE
  diff_in_means <- subset(FoD, D_i==1)$diff_mean-subset(FoD, D_i==0)$diff_mean
  se_diff_in_means <- sqrt(subset(FoD, D_i==1)$s2*(subset(FoD, D_i==1)$n-1)/(subset(FoD, D_i==1)$n)^2+subset(FoD, D_i==0)$s2*(subset(FoD, D_i==0)$n-1)/(subset(FoD, D_i==0)$n)^2) # var() uses denominator (n_g-1); i'd like n_g for comparison with HC0
  
  cov_diff_in_means <- c(pl_hldr, true_tau >= diff_in_means - 1.96*se_diff_in_means & true_tau <= diff_in_means + 1.96*se_diff_in_means)
  diff_in_means <- c(pl_hldr, diff_in_means)
  se_diff_in_means <- c(pl_hldr, se_diff_in_means)
  names(se_diff_in_means) <- c(1979:2016, 2018)
  names(diff_in_means) <- c(1979:2016, 2018)
  names(cov_diff_in_means) <- c(1979:2016, 2018)
  
  # run two-group, multiple pre-period event study
  time_start <- Sys.time()
  lm_ES <- lm(Y_i ~ D_i+T_f+D_i:T_f, data=pnl_data_noSDID)
  
  # get heteroskedasticity-robust (Eicker–Huber–White) SE clustered at the person level for panel data (independent across i and states, with serial correlation across time for each i)
  SE <- vcovCL(lm_ES, cluster = ~person_id, type = "HC0")
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
    c(as.list(tau_pt), list(type = "est", DGP = 4, niter=iter)),
    c(as.list(se_pt), list(type = "se", DGP = 4, niter=iter)),
    c(as.list(coverage_pt), list(type = "coverage", DGP = 4, niter=iter)),
    c(as.list(time_elapsed_pt), list(type = "time", DGP = 4, niter=iter)),
    c(as.list(diff_in_means), list(type = "est_dim", DGP = 4, niter=iter)),
    c(as.list(se_diff_in_means), list(type = "se_dim", DGP = 4, niter=iter)),
    c(as.list(cov_diff_in_means), list(type = "coverage_dim", DGP = 4, niter=iter))
  )
  ## END: PT-DID ------------------------------------
  
  
  ## START: SDID ------------------------------------
  Y <- L + rmvnorm(50, sigma = Sigma, method = "chol")
  Y <- Y[order(trt_idx), ]
  
  time_start <- Sys.time()
  method_sdid <- synthdid_estimate(Y, 
                                   factor_model_1tr$N0,
                                   factor_model_1tr$T0)
  tau_sdid <- summary(method_sdid)$estimate
  se_sdid <- sqrt(vcov(method_sdid, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sdid <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sdid <- as.numeric(true_tau >= tau_sdid - 1.96*se_sdid & true_tau <= tau_sdid + 1.96*se_sdid)
  
  # collect results
  result_sdid <- bind_rows(
    result_sdid,
    list("2018" = tau_sdid, type = "est", DGP = 4, niter=iter),
    list("2018" = se_sdid, type = "se", DGP = 4, niter=iter),
    list("2018" = coverage_sdid, type = "coverage", DGP = 4, niter=iter),
    list("2018" = time_elapsed_sdid, type = "time", DGP = 4, niter=iter)
  )
  ## END: SDID --------------------------------------
  
  ## START: SC ------------------------------------
  time_start <- Sys.time()
  
  method_sc <- sc_estimate(Y, 
                           factor_model_1tr$N0,
                           factor_model_1tr$T0)
  tau_sc <- summary(method_sc)$estimate
  se_sc <- sqrt(vcov(method_sc, method='placebo'))
  
  time_end <- Sys.time()
  time_elapsed_sc <- as.numeric(difftime(time_end, time_start, units = "secs"))
  
  coverage_sc <- as.numeric(true_tau >= tau_sc - 1.96*se_sc & true_tau <= tau_sc + 1.96*se_sc)
  
  # collect results
  result_sc <- bind_rows(
    result_sc,
    list("2018" = tau_sc, type = "est", DGP = 4, niter=iter),
    list("2018" = se_sc, type = "se", DGP = 4, niter=iter),
    list("2018" = coverage_sc, type = "coverage", DGP = 4, niter=iter),
    list("2018" = time_elapsed_sc, type = "time", DGP = 4, niter=iter)
  )
  ## END: SC --------------------------------------
  
  
  ## START: SPT ------------------------------------
  time_start <- Sys.time()
  n <- length(pnl_data_noSDID$Y_i)
  
  method_spt <- tau_set(outcome=pnl_data_noSDID$Y_i,
                        unit_id=pnl_data_noSDID$G_i,
                        time_id=pnl_data_noSDID$T_i,
                        trt_unit_id=which(trt_idx==1)+50,
                        person_id=pnl_data_noSDID$person_id,
                        cand_tau_size=cand_tau_size,
                        num_bstp_rep=num_bstp_rep,
                        cand_tau=cand_tau_noSDID)
  
  time_end <- Sys.time()
  time_elapsed_spt <- as.numeric(difftime(time_end, time_start, units = "secs"))
  cs_tau <- method_spt$cand_tau[method_spt$rej==0]
  
  # for future iterations, candidate values of tau for the next iteration are based on enlarging current iteration's cs_tau, if nonempty
  if (length(cs_tau)>0){
    # double the distance from the midpoint (symmetric expansion)
    a <- min(cs_tau)
    b <- max(cs_tau)
    mid <- (a+b)/2
    len <- b-a
    a <- mid-len
    b <- mid+len
    cand_tau_noSDID <- a+(b-a)*
      qbeta(seq(0, 1, length.out = cand_tau_size), shape1 = 2, shape2 = 2) 
  }
  
  # collect results
  result_spt <- bind_rows(
    result_spt,
    list("2018" = ifelse(length(cs_tau)>0, min(cs_tau), -99), type = "cs_lb", DGP = 4, niter=iter),
    list("2018" = ifelse(length(cs_tau)>0, max(cs_tau), -99), type = "cs_ub", DGP = 4, niter=iter),
    list("2018" = ifelse(length(cs_tau)>0, as.numeric(true_tau <= max(cs_tau) & true_tau >= min(cs_tau)), 0), type = "coverage", DGP = 4, niter=iter),
    list("2018" = time_elapsed_spt, type = "time", DGP = 4, niter=iter)
  )
  frac_rej_noSDID <- append(frac_rej_noSDID, 1-length(cs_tau)/cand_tau_size)
  ## END: SPT ------------------------------------
  ### END: DGP-4 treated unit is one state, many weights are valid pre but only a sparse weight is valid post ----
  
  end <- Sys.time()
  
  print(paste0("DGP-3: Many valid weights pre AND post ------------"))
  print(paste0("PT-DID avg bias: ", mean(abs(true_tau-subset(result_pt, type=="est" & DGP==3)$`2018`))))
  print(paste0("PT-DID avg bias (DiM): ", mean(abs(true_tau-subset(result_pt, type=="est_dim" & DGP==3)$`2018`))))
  print(paste0("PT-DID avg CI length: ", 1.96*2*mean(subset(result_pt, type=="se" & DGP==3)$`2018`)))
  print(paste0("PT-DID avg CI length (DiM): ", 1.96*2*mean(subset(result_pt, type=="se_dim" & DGP==3)$`2018`)))
  print(paste0("PT-DID coverage rate: ", mean(subset(result_pt, type=="coverage" & DGP==3)$`2018`)))
  print(paste0("PT-DID coverage rate (DiM): ", mean(subset(result_pt, type=="coverage_dim" & DGP==3)$`2018`)))
  print(paste0("PT-DID avg time: ", mean(subset(result_pt, type=="time" & DGP==3)$`2018`)))
  print(paste0("-"))
  print(paste0("SDID avg bias: ", mean(abs(true_tau-subset(result_sdid, type=="est" & DGP==3)$`2018`))))
  print(paste0("SDID avg CI length: ", 1.96*2*mean(subset(result_sdid, type=="se" & DGP==3)$`2018`)))
  print(paste0("SDID coverage rate: ", mean(subset(result_sdid, type=="coverage" & DGP==3)$`2018`)))
  print(paste0("SDID avg time: ", mean(subset(result_sdid, type=="time" & DGP==3)$`2018`)))
  print(paste0("-"))
  print(paste0("SC avg bias: ", mean(abs(true_tau-subset(result_sc, type=="est" & DGP==3)$`2018`))))
  print(paste0("SC avg CI length: ", 1.96*2*mean(subset(result_sc, type=="se" & DGP==3)$`2018`)))
  print(paste0("SC coverage rate: ", mean(subset(result_sc, type=="coverage" & DGP==3)$`2018`)))
  print(paste0("SC avg time: ", mean(subset(result_sc, type=="time" & DGP==3)$`2018`)))
  print(paste0("-"))
  print(paste0("SPT avg CS length: ", mean(subset(result_spt, type=="cs_ub" & DGP==3)$`2018`-subset(result_spt, type=="cs_lb" & DGP==3)$`2018`)))
  print(paste0("SPT coverage rate: ", mean(subset(result_spt, type=="coverage" & DGP==3)$`2018`)))
  print(paste0("SPT avg time: ", mean(subset(result_spt, type=="time" & DGP==3)$`2018`)))
  print(paste0("frac cand_tau rej: ", mean(frac_rej_SDID)))
  print(paste0("DGP-4: Many valid weights pre but NOT post ------------"))
  print(paste0("PT-DID avg bias: ", mean(abs(true_tau-subset(result_pt, type=="est" & DGP==4)$`2018`))))
  print(paste0("PT-DID avg bias (DiM): ", mean(abs(true_tau-subset(result_pt, type=="est_dim" & DGP==4)$`2018`))))
  print(paste0("PT-DID avg CI length: ", 1.96*2*mean(subset(result_pt, type=="se" & DGP==4)$`2018`)))
  print(paste0("PT-DID avg CI length (DiM): ", 1.96*2*mean(subset(result_pt, type=="se_dim" & DGP==4)$`2018`)))
  print(paste0("PT-DID coverage rate: ", mean(subset(result_pt, type=="coverage" & DGP==4)$`2018`)))
  print(paste0("PT-DID coverage rate: ", mean(subset(result_pt, type=="coverage_dim" & DGP==4)$`2018`)))
  print(paste0("PT-DID avg time: ", mean(subset(result_pt, type=="time" & DGP==4)$`2018`)))
  print(paste0("-"))
  print(paste0("SDID avg bias: ", mean(abs(true_tau-subset(result_sdid, type=="est" & DGP==4)$`2018`))))
  print(paste0("SDID avg CI length: ", 1.96*2*mean(subset(result_sdid, type=="se" & DGP==4)$`2018`)))
  print(paste0("SDID coverage rate: ", mean(subset(result_sdid, type=="coverage" & DGP==4)$`2018`)))
  print(paste0("SDID avg time: ", mean(subset(result_sdid, type=="time" & DGP==4)$`2018`)))
  print(paste0("-"))
  print(paste0("SC avg bias: ", mean(abs(true_tau-subset(result_sc, type=="est" & DGP==4)$`2018`))))
  print(paste0("SC avg CI length: ", 1.96*2*mean(subset(result_sc, type=="se" & DGP==4)$`2018`)))
  print(paste0("SC coverage rate: ", mean(subset(result_sc, type=="coverage" & DGP==4)$`2018`)))
  print(paste0("SC avg time: ", mean(subset(result_sc, type=="time" & DGP==4)$`2018`)))
  print(paste0("-"))
  print(paste0("SPT avg CS length: ", mean(subset(result_spt, type=="cs_ub" & DGP==4)$`2018`-subset(result_spt, type=="cs_lb" & DGP==4)$`2018`)))
  print(paste0("SPT coverage rate: ", mean(subset(result_spt, type=="coverage" & DGP==4)$`2018`)))
  print(paste0("SPT avg time: ", mean(subset(result_spt, type=="time" & DGP==4)$`2018`)))
  print(paste0("frac cand_tau rej: ", mean(frac_rej_noSDID)))
  
  print(end-start)
  
  write.csv(result_spt, "../result/result_spt_DGP34.csv", row.names=FALSE)
  write.csv(result_pt, "../result/result_pt_DGP34.csv", row.names=FALSE)
  write.csv(result_sdid, "../result/result_sdid_DGP34.csv", row.names=FALSE)
  write.csv(result_sc, "../result/result_sc_DGP34.csv", row.names=FALSE)
}
