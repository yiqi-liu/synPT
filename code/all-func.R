# generate fake repeated cross-sectional sample for control units
# DGP-1: PT==TRUE  & SDID==FALSE
# DGP-2: PT==FALSE & SDID==FALSE
# DGP-3: PT==FALSE & SDID==TRUE
gen_micro_data_control <- function(k, t=NULL, n, 
                                   PT=TRUE, SDID=FALSE, Sigma){
  if (PT & SDID) stop(paste0("Only one DGP specification is allowed."))
  
  if (!SDID){
    if (!PT & t==2018){ # if DGP-2, we generate post-treatment outcomes differently
      m <- ctrl_mean_post_noPT[rownames(ctrl_mean_post_noPT)==as.character(k)]
    } else {
      m <- subset(cell_stats, state_num==k & year==t)$wage 
    }
    s <- subset(cell_stats, state_num==k & year==t)$sd
    
    tibble(G_i=k, T_i=t, D_i=0,
           Y_i=rnorm(n, m, s))
  } else { # create panel for DGP-3
    Y_it <- rmvnorm(n, mean = L[k, ], sigma = n*Sigma, method = "chol") # n by 40
    colnames(Y_it) <- 1979:2018
    
    dta <- as_tibble(Y_it) %>%
      mutate(G_i=k, D_i=0, person_id=1:n) %>%
      pivot_longer(cols=`1979`:`2018`, names_to="T_i", values_to="Y_i") %>%
      mutate(T_i = as.numeric(T_i))
    
    return(dta)
  }
}

# generate population means for the treated unit; only used for DGP-1 and DGP-2
gen_trt_unit <- function(trt_idx, pi_0=NULL, PT=TRUE){
  # for the DGP satisfying PT, the period-T counterfactual mean outcome of the simulated treated unit is the average of the period-T means of the selected n_1 treated state(s) from calling 'simulate_dgp'
  trt_counterfactual <- mean(subset(cell_stats, state_num %in% which(trt_idx==1) & year==2018)$wage)
  
  if (is.null(pi_0)){ # weights used to generate trends of the treated unit; if not specified, take the weighted average of control trends, weighted by pi_0k, the control frequencies in the last period
    pi_0 <- subset(cell_stats, year==2018 & state_num %in% which(trt_idx==0))$n_kt/sum(subset(cell_stats, year==2018 & state_num %in% which(trt_idx==0))$n_kt)  
  } 
  
  # do the same for the selected treated states
  pi_1 <- subset(cell_stats, year==2018 & state_num %in% which(trt_idx==1))$n_kt/sum(subset(cell_stats, year==2018 & state_num %in% which(trt_idx==1))$n_kt)
  
  # calculate trends of the treated using pi
  trt_trends <- c()
  for (yr in 1980:2018){
    trt_trends <- append(trt_trends,
                         sum(pi_0*(subset(cell_stats, year==yr & state_num %in% which(trt_idx==0))$wage-subset(cell_stats, year==yr-1 & state_num %in% which(trt_idx==0))$wage))
    )
  }
  
  if (!PT){ # generate nonparallel post-trend in the noPT DGP 
    trt_trends[length(trt_trends)] <- sum(as.numeric(a_post_alt)*as.numeric(w_alt))
  }
  
  # generate the untreated potential means for the treated unit
  trt_mean <- rep(0, 40)
  trt_mean[40] <- trt_counterfactual
  for (yr in 39:1){
    trt_mean[yr] <- trt_mean[yr+1]-trt_trends[yr]
  }
  
  return(list(trt_mean=trt_mean,
              trt_trends=trt_trends,
              pi_0=pi_0, pi_1=pi_1))
}

# generate repeated cross-sectional sample for the treated unit
gen_micro_data_trt <- function(trt_info=NULL, k, t=NULL, n, true_tau, 
                               SDID=FALSE, Sigma){
  if (!SDID){
    if (t==2018){ # post period, will add true_tau
      m <- subset(cell_stats, state_num==k & year==t)$wage
      s <- subset(cell_stats, state_num==k & year==t)$sd
      
      tibble(G_i=50+k, T_i=t, D_i=1,  # assign temporary group id=50+k
             Y_i=rnorm(n, m+true_tau, s))
    } else{
      mu_avg <- trt_info$trt_mean[t-1979+1]
      m <- sum(trt_idx)*mu_avg*n/n_1
      s <- sd(subset(rcs_CPS, state_num==k & year==t)$outcome)  # pooled sd
      
      tibble(G_i=50+k, T_i=t, D_i=1, # assign temporary group id=50+k
             Y_i=rnorm(n, m, s))
    }
  } else{ # create panel for DGP-3
    L_mod <- L
    L_mod[k, 40] <- L_mod[k, 40] + true_tau # add tau to last year
    
    Y_it <- rmvnorm(n, mean = L_mod[k, ], sigma = n*Sigma, method = "chol") # n by 40
    
    colnames(Y_it) <- 1979:2018
    
    dta <- as_tibble(Y_it) %>%
      mutate(G_i=50+k, D_i=1, person_id=1:n) %>% # assign temporary group id=50+k
      pivot_longer(cols=`1979`:`2018`, names_to="T_i", values_to="Y_i") %>%
      mutate(T_i=as.numeric(T_i))
    
    return(dta)
  }
}


# function that constructs a confidence set for the identified set of treatment effect on the treated unit in the post-treatment period
tau_set <- function(outcome,  # 1 by n vector of outcomes
                    unit_id,  # 1 by n vector of unit ids; should be numeric
                    time_id,  # 1 by n vector of year ids; should be numeric; time periods should already be sorted in increasing order with the difference between each consecutive period set to 1; the largest number is taken to be the post-treatment period
                    trt_unit_id, # the id of the treated unit; should be numeric
                    person_id=NULL, # used for panel data; if NULL, assume repeated cross sections
                    cand_tau=NULL, # candidate values for tau to test invert
                    cand_tau_size=100, # number of candidate values for tau
                    num_bstp_rep=NULL, # number of bootstrap replications; if NULL, then just return the test statistic without bootstrapped critical values
                    s_n=NULL, # step size for numerical approximation of the directional derivative
                    alpha=0.05, # significance level; not used if num_bstp_rep=NULL
                    infntsml_adj=1e-6, # infinitesimal adjustment factor; not used if num_bstp_rep=NULL
                    tol=NULL # tolerance to abs(A_pre w - b_pre) <= tol
){
  K0 <- length(unique(unit_id))-1 # number of control units
  T0 <- length(unique(time_id))-1 # number of pre-periods
  post_T <- max(time_id)  # post-treatment period
  n <- length(outcome)    # effective sample size

  # get kt cell means
  cell_means <- data.frame(outcome=outcome,
                           unit_id=unit_id,
                           time_id=time_id) %>%
    group_by(unit_id, time_id) %>%
    summarise(mu=mean(outcome, na.rm = TRUE), .groups = "drop")
  
  # create treatment indicator
  cell_means <- cell_means %>% 
    mutate(W = unit_id==trt_unit_id & time_id==post_T)
  
  # get a K0+1 by T0+1 outcome matrix using 'panel.matrices' from synthdid, where the outcomes for the treated unit is sorted in the last row
  outcome_mat <- synthdid::panel.matrices(data.frame(cell_means))$Y
  
  # get a K0+1 by T0 trend matrix
  trend_mat <- outcome_mat[, 2:(T0+1), drop=FALSE]-outcome_mat[, 1:T0, drop=FALSE]
  
  # construct A_pre (T_0-1 by K_0)
  A_pre <- trend_mat[1:K0, 1:(T0-1), drop=FALSE] # K0 by (T0-1)
  A_pre <- t(A_pre[order(as.numeric(row.names(A_pre))), , drop=FALSE]) # (T0-1) by K0; control indices are sorted in increasing order
  
  # construct b_pre (T_0-1 by 1)
  b_pre <- matrix(trend_mat[K0+1, 1:(T0-1), drop=FALSE])
  
  # construct a_post (K0 by 1)
  a_post <- matrix(trend_mat[1:K0, T0, drop=FALSE])
  
  # construct the observed difference in period-T and period-T0 mean outcome of the treated unit
  d_post <- as.numeric(trend_mat[K0+1, T0]) # scalar            
  
  # Rewrite (Aw-b)'(Aw-b)+(aw-d+tau) = w'(A'A+aa')w + w'(-2A'b-2ad+2a*tau) + (const independent of w) = w'Pw/2 + (q_base+2a*tau)'w + const for P = 2(A'A+aa'); q_base = 2(-A'b-ad); let q(tau) = q_base + 2a*tau
  P  <- 2*(crossprod(A_pre)+tcrossprod(a_post)) # K0 x K0
  q_base <- 2*(-crossprod(A_pre, b_pre)-as.numeric(a_post)*d_post)  # K0 x 1 numeric
  
  # simplex constraint on w: w >= 0, sum(w) == 1
  # in terms of l <= Cw <= u
  simplex_constr <- rbind(Diagonal(K0), 
                          Matrix(1, nrow = 1, ncol = K0, 
                                 sparse = TRUE))
  l <- c(rep(0, K0), 1); u <- c(rep(Inf, K0), 1)
  
  # build OSQP once that does not depend on tau
  osqp_solver <- osqp(
    P = as(P, "dgCMatrix"),
    q = as.numeric(q_base),
    A = as(simplex_constr, "dgCMatrix"),
    l = l, u = u,
    osqpSettings(verbose = FALSE, polish=TRUE,
                 adaptive_rho_interval = 50 # for reproducibility; see https://github.com/osqp/osqp-r/issues/19#issuecomment-636954982
                 )
  )
  
  # create candidate values for tau if not specified
  if (is.null(cand_tau)){
    cand_tau <- gen_cand_tau(A_pre, b_pre, a_post, d_post, 
                             size=cand_tau_size, 
                             tol=tol)
  }
  
  # solve for each tau by updating q(tau) only
  test_stat <- numeric(length(cand_tau))
  w_prev <- rep(1/K0, K0) # initial guess is equal weight
  for (j in seq_along(cand_tau)) {
    tau <- cand_tau[j] # current tau value
    # update q(tau)
    q_tau <- as.numeric(q_base+2*as.numeric(a_post)*tau) 
    osqp_solver$Update(q = q_tau)
    osqp_solver$WarmStart(x = w_prev)
    # solve
    sol <- osqp_solver$Solve()
    if (!(sol$info$status_val %in% c(1L, 2L))){
      stop(paste0("OSQP status: ", sol$info$status))
    }
    
    # get the corresponding optimal w for current tau
    w_opt <- sol$x
    # update initial solution
    w_prev <- w_opt
    # evaluate objective value (Aw-b)'(Aw-b)+(aw-d+tau) at current optimal w; directly extracting optimal value from sol would be wrong because the QP objective ignores the part of (Aw-b)'(Aw-b)+(aw-d+tau) that doesn't depend on w
    r_pre <- A_pre%*%w_opt-b_pre
    r_post <- as.numeric(crossprod(a_post, w_opt)-(d_post - tau))
    test_stat[j] <- sqrt(n)*(sum(r_pre^2) + r_post^2) # compute test statistic
  }
  
  # initialize vector to collect bootstrapped cv
  btstrp_cv <- c()
  if (!is.null(num_bstp_rep)) {
    # step size for approximating the directional derivative of min_{convex \omega}(⋅), as per Eq. (25) of Fang & Santos (2010)
    if (is.null(s_n)) s_n <- n^(-1/2+0.01)
    # initialize matrix to store bootstrapped test statistic
    BS_test_stat <- matrix(nrow=length(cand_tau), ncol=num_bstp_rep)
    
    # mapping from cand_tau to indices
    tau_to_idx <- setNames(seq_along(cand_tau), as.character(cand_tau))
    
    # bootstrap once for all param
    ### ------- BOOTSTRAP STARTS -------
    for (bs in seq_len(num_bstp_rep)) {
      # bootstrapped kt cell means
      if (is.null(person_id)){  # without person_id, assume data is repeated cross-section
        # draw n exponential(1) weights
        bs_weight <- rexp(n, rate = 1)
        # reweight the cell means
        bs_cell_means <- data.table(
          unit_id=unit_id,
          time_id=time_id,
          outcome=outcome, 
          weight=bs_weight
        )[, .(bs_mean=sum(weight*outcome)/sum(weight)), 
          by = .(unit_id, time_id)]  %>% 
          mutate(W = unit_id==trt_unit_id & time_id==post_T)
        
      } else{ # if it's panel data
        pid_key <- interaction(unit_id, person_id, drop = TRUE) # unit-person key; length=n
        pid_u <- levels(pid_key) # unique persons i across all aggregate units; length=n/40
        bs_w_u <- rexp(length(pid_u), rate = 1) # each person i gets the same weight used for all periods
        bs_weight <- bs_w_u[pid_key] # fill in length=n vector of weights
        
        bs_cell_means <- data.table(
          unit_id = unit_id,
          time_id = time_id,
          outcome = outcome,
          weight  = bs_weight
        )[ , .(bs_mean = sum(weight * outcome) / sum(weight)),
           by = .(unit_id, time_id) ] %>% 
          mutate(W = unit_id==trt_unit_id & time_id==post_T)
      }
      
      # follows the same steps as before to construct bootstrapped versions of A_pre, b_pre, a_post, and d_post
      bs_outcome_mat <- synthdid::panel.matrices(data.frame(bs_cell_means))$Y
      bs_trend_mat <- bs_outcome_mat[, 2:(T0+1), drop=FALSE] - bs_outcome_mat[, 1:T0, drop=FALSE]
      
      bs_A_pre <- bs_trend_mat[1:K0, 1:(T0-1), drop=FALSE]
      bs_A_pre <- t(bs_A_pre[order(as.numeric(rownames(bs_A_pre))), , drop=FALSE])
      bs_b_pre <- matrix(bs_trend_mat[K0+1, 1:(T0-1), drop=TRUE])
      bs_a_post <- matrix(bs_trend_mat[1:K0, T0, drop=TRUE])
      bs_d_post <- as.numeric(bs_trend_mat[K0+1, T0])
      
      # compute bootstrapped matrices: original + s_n*sqrt(n)*direction
      eval_A_pre <- A_pre+s_n*sqrt(n)*(bs_A_pre-A_pre)
      eval_b_pre <- b_pre+s_n*sqrt(n)*(bs_b_pre-b_pre)
      eval_a_post <- a_post+s_n*sqrt(n)*(bs_a_post-a_post)
      eval_d_post <- d_post+s_n*sqrt(n)*(bs_d_post-d_post)
      
      # build eval_P and eval_q_base once per bootstrap, reuse across tau
      eval_P <- 2*(crossprod(eval_A_pre)+tcrossprod(eval_a_post))
      eval_q_base <- 2*(-crossprod(eval_A_pre, eval_b_pre)-as.numeric(eval_a_post)*eval_d_post)
      
      bs_osqp <- osqp(
        P = as(eval_P, "dgCMatrix"),
        q = as.numeric(eval_q_base),
        A = as(simplex_constr, "dgCMatrix"),
        l = l, u = u,
        osqpSettings(verbose = FALSE,
                     adaptive_rho_interval = 50 # for reproducibility; see https://github.com/osqp/osqp-r/issues/19#issuecomment-636954982
                     )
      )
      
      # initialize vector to store bootstrapped test statistic
      bs_vec <- numeric(length(cand_tau))
      bs_w_prev <- w_prev # initial value guess
      
      # loop over candidate tau
      for (j in seq_along(cand_tau)) {
        tau <- cand_tau[j]
        eval_q_tau <- as.numeric(eval_q_base+2*as.numeric(eval_a_post)*tau) # get q(tau) for perturbed estimators
        bs_osqp$Update(q=eval_q_tau)
        bs_osqp$WarmStart(x=bs_w_prev)
        
        solb <- bs_osqp$Solve()
        if (!(solb$info$status_val %in% c(1L, 2L))){
          stop(paste0("OSQP (bootstrap) status: ", solb$info$status))
        } 
        
        bs_w_opt <- solb$x
        bs_w_prev <- bs_w_opt
        
        # evaluate objective value (Aw-b)'(Aw-b)+(aw-d+tau) at optimal bs_w_opt using perturbed estimators:
        bs_r_pre <- eval_A_pre%*%bs_w_opt-eval_b_pre
        bs_r_post <- as.numeric(crossprod(eval_a_post, bs_w_opt)-(eval_d_post-tau))
        eval_phi <- sum(bs_r_pre^2)+bs_r_post^2
        
        # original objective value
        orig_phi <- test_stat[j]/sqrt(n)
        
        # get bootstrapped test statistic
        bs_vec[j] <- (eval_phi-orig_phi)/s_n
      }
      
      BS_test_stat[, bs] <- bs_vec
    }
    
    # bootstrap critical values at the (1-alpha) level
    btstrp_cv <- apply(
      BS_test_stat, 1L,
      function(x) quantile(x, 1-alpha+infntsml_adj, na.rm=TRUE)+infntsml_adj
    )
  }
  
  list(
    cand_tau = cand_tau,
    test_stat = test_stat,
    btstrp_cv = btstrp_cv, # will return empty c() if num_bstp_rep=NULL
    rej = as.numeric(test_stat > btstrp_cv)
  )
}



# function that creates candidate values for tau by running an initial LP
gen_cand_tau <- function(A_pre, # (T_0-1 by K_0) control pre-trends matrix
                         b_pre, # (T_0-1 by 1) treated pre-trends matrix
                         a_post, # (K0+1 by 1) control post-trend vector
                         d_post, # scalar, observed post-trend of the treated
                         size=100, # output number of candidate values
                         tol=NULL # tolerance to abs(A_pre w - b_pre) <= tol
){
  b_pre <- as.numeric(b_pre)
  a_post <- as.numeric(a_post)
  
  K0 <- ncol(A_pre)
  T0 <- nrow(A_pre)+1
  
  # build LP: min a_post'w subject to w: w >= 0; 1'w = 1;  A_pre w - b_pre <= tol;  -A_pre w + b_pre <= tol
  # simplex constraint:
  constr_simplex <- rbind(Diagonal(K0), 
                          Matrix(1, nrow = 1, ncol = K0, sparse = TRUE))
  l_simplex <- c(rep(0, K0), 1)
  u_simplex <- c(rep(Inf, K0), 1)
  
  if (is.null(tol)){
    # set tol to be min_{w convex} (1/2)||A_pre w - b_pre||^2 = (1/2) w'(A_pre'A_pre)w-b_pre'A_pre w + const 
    QP_tol <- osqp(
      P = as(crossprod(A_pre), "dgCMatrix"), # A_pre'A_pre
      q = as.numeric(-crossprod(A_pre, b_pre)), # -A_pre'b_pre
      A = constr_simplex,
      l = l_simplex, u = u_simplex,
      osqpSettings(verbose = FALSE, polish = TRUE,
                   adaptive_rho_interval = 50 # for reproducibility; see https://github.com/osqp/osqp-r/issues/19#issuecomment-636954982
                   )
    )
    
    res_tol <- QP_tol$Solve()
    w_tol <- res_tol$x
    tol <- norm(A_pre %*% w_tol-b_pre, type = "2")
  }

  # pre-trt constraints
  # A_pre w - b_pre <= tol
  pre_constr <- Matrix(A_pre, sparse = TRUE)
  l_pre_constr <- rep(-Inf, T0-1)
  u_pre_constr <- b_pre + tol
  
  # -A_pre w + b_pre <= tol
  neg_pre_constr <- -pre_constr
  l_neg_pre <- rep(-Inf, T0-1)
  u_neg_pre <- -b_pre + tol
  
  # combine
  constr <- rbind(constr_simplex, pre_constr, neg_pre_constr)
  l <- c(l_simplex, l_pre_constr, l_neg_pre)
  u <- c(u_simplex, u_pre_constr, u_neg_pre)
  
  # min (1/2) w'Pw + q'w with P = 0
  LP <- osqp(P = Diagonal(K0), # LP has no quadratic component
             q = a_post, # min a_post'w
             A = as(constr, "dgCMatrix"),
             l = as.numeric(l),
             u = as.numeric(u),
             osqpSettings(verbose = FALSE, polish = TRUE, 
                          max_iter = 200000,
                          adaptive_rho_interval = 50 # for reproducibility; see https://github.com/osqp/osqp-r/issues/19#issuecomment-636954982
                          )
             )
  
  # min a_post'w
  res_min <- LP$Solve()
  if (!res_min$info$status_val %in% c(1L, 2L)) stop(sprintf("OSQP(min) status: %s", res_min$info$status))
  val_min <- as.numeric(res_min$info$obj_val) 
  
  # max a_post'w
  LP$Update(q=-a_post)
  LP$WarmStart(x=res_min$x)
  res_max <- LP$Solve()
  if (!res_max$info$status_val %in% c(1L, 2L)) stop(sprintf("OSQP(max) status: %s", res_max$info$status))
  val_max <- -as.numeric(res_max$info$obj_val)
  
  tau_lb <- d_post-val_max
  tau_ub <- d_post-val_min
  
  # enlarge the distance from the midpoint
  mid <- (tau_lb+tau_ub)/2
  len <- tau_ub-tau_lb
  tau_lb <- mid-2*len
  tau_ub <- mid+2*len
  
  # sample more points in the interior of [tau_lb, tau_ub]
  cand_tau <- tau_lb+(tau_ub-tau_lb)*
    qbeta(seq(0, 1, length.out = size), shape1 = 2, shape2 = 2)
  
  return(cand_tau)
}


## function that finds a different convex weight that generates parallel pre-trends but not parallel post-trends:
get_alt_omega <- function(A_pre, b_pre, a_post, pi_0, 
                          delta_target, # post-trend violation
                          eps_reg = 1e-8, verbose = FALSE) {

  b_pre  <- as.numeric(b_pre)
  a_post <- as.numeric(a_post)
  pi_0    <- as.numeric(pi_0)
  K0 <- ncol(A_pre)
  
  # parallel post-trends are generated by weighting a_post by pi_0
  rhs_base <- as.numeric(crossprod(a_post, pi_0))
  
  # base constraints: A_pre w = b_pre; sum(w) = 1; w >= 0
  constr_base <- rbind(
    Matrix(A_pre, sparse = TRUE),
    Matrix(1, nrow = 1, ncol = K0, sparse = TRUE),
    Diagonal(K0)
  )
  l_base <- c(b_pre, 1, rep(0, K0))
  u_base <- c(b_pre, 1, rep(Inf, K0))
  
  # small reg to make QP strictly convex; objective is close to 0 so just checking feasibility
  P <- as(eps_reg * Diagonal(K0), "dgCMatrix")
  q <- rep(0, K0)
  
  solve_dir <- function(direction = c("up", "down")) {
    direction <- match.arg(direction)
    
    if (direction == "up") {
      # add constraint: a_post' w >= rhs_base + delta_target
      constr <- rbind(constr_base, Matrix(rbind(a_post), sparse = TRUE))
      l <- c(l_base, rhs_base+delta_target)
      u <- c(u_base, Inf)
    } else {
      # add constraint: a_post' w <= rhs_base - delta_target
      constr <- rbind(constr_base, Matrix(rbind(a_post), sparse = TRUE))
      l <- c(l_base, -Inf)
      u <- c(u_base, rhs_base-delta_target)
    }
    
    model <- osqp(
      P = P,
      q = as.numeric(q),
      A = as(constr, "dgCMatrix"),
      l = as.numeric(l),
      u = as.numeric(u),
      osqpSettings(verbose = verbose, polish = TRUE,
                   eps_abs = 1e-8, eps_rel = 1e-8, max_iter = 200000,
                   adaptive_rho_interval = 50 # for reproducibility; see https://github.com/osqp/osqp-r/issues/19#issuecomment-636954982
                   )
    )
    res <- model$Solve()
    if (!(res$info$status_val %in% c(1L, 2L))) return(NULL)
    
    w <- res$x
    delta <- abs(as.numeric(crossprod(a_post, w-pi_0)))
    
    if (round(delta-delta_target, 5) >= 0) return(list(w = w, delta = delta))
    return(NULL)
  }
  
  # try "up" first
  res_up <- solve_dir("up")
  if (!is.null(res_up)) {
    return(list("w_alt" = res_up$w, "violation"=res_up$delta, "a_post_alt"=a_post))
  }
  
  # try "down" if "up" fails
  res_down <- solve_dir("down")
  if (!is.null(res_down)) {
    return(list("w_alt" = res_down$w, "violation"=res_down$delta, "a_post_alt"=a_post))
  }

  return(NULL)
}




