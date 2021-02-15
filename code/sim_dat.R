
sample_generator = function(n, C_delta, treatment_type, control_type = "bell",
                            m = NA, level = NA, paired = FALSE,  eps = NA){
  dat = list()
  A = rbinom(n = n, size = 1, prob = 1/2)
  
  if (any(grepl("subgroup", treatment_type))) {
    X = cbind(sample(1:level, n, replace = TRUE),
              sample(1:2, n, replace = TRUE))
  } else {
    X = cbind(c(rep(0, n/2), rep(1, n/2)), #a1
              c(rep(0, m), rep(1, n/2 - m), rep(0, n/2-m), rep(1, m)), #a2
              rnorm(n = n, sd = 1)) #a3 
  }
  
  if (treatment_type == "sparse_pos_bias") {
    Delta = function(X) {5*(X[,3])^3*(X[,3] > 1) - X[,1]/2}
  } else if (treatment_type == "linear") { 
    Delta = function(X) {2*(X[,1]*X[,2] + X[,3])}
  } else if (treatment_type == "sparse_oneside") {
    Delta = function(X) {5*(X[,3])^3*(X[,3] > 1)}
  } else if (treatment_type == "sparse_twoside") {
    Delta = function(X) {5*(X[,3])^3*(abs(X[,3]) > 1)}
  } else if (treatment_type == "subgroup_even") {
    # if (paired) {
    #   Delta = function(X) {(X[,1] %% 2 == 0)}
    # } else {
    #   Delta = function(X) {0.2*(X[,1] %% 2 == 0)}
    # }
    Delta = function(X) {X[,1] %% 2 == 0}
  } else if (treatment_type == "subgroup_smooth") {
    Delta = function(X) {X[,1] <= 20} 
  } else if (treatment_type == "subgroup_sparse") {
    Delta = function(X) {X[,1] %% 4 == 0}
  } 

  if (control_type == "bell"){
    f = function(X) {5*rowSums(X)}
  } else if (control_type == "skewed") {
    f = function(X) {2*((X[,3] < -2)*exp(-X[,3])^2)}
  } else if (control_type == "zero") {
    f = function(X) {0}
  } else if (control_type == "nonlinear") {
    f = function(X) {5*abs(rowSums(X))}
  } 

  if(!paired) {
    Y = A*C_delta*Delta(X) + f(X) + rnorm(n)
    dat$prob = rep(1/2, n)
    dat$observe = list(Y = Y, A = A, X = X)
    dat$nonnull_ind = C_delta*Delta(X) != 0
    dat$nonnull_pos = C_delta*Delta(X) > 0
  } else {
    A1 = A; A2 = 1 - A
    X1 = X2 = X;
    if(!is.na(eps)) {
      X2[,3] = X[,3] + runif(n, min = 0, max = 2*eps)
      eps = min(1, eps)
      change_one = rbinom(n, 1, eps); X2[,1] = (1 - X[,1])*change_one + X[,1]*(1 - change_one)
      change_two = rbinom(n, 1, eps); X2[,2] = (1 - X[,2])*change_two + X[,2]*(1 - change_two)
    }
    
    Y1 = A1*C_delta*Delta(X1) + f(X1) + rnorm(n)
    Y2 = A2*C_delta*Delta(X2) + f(X2) + rnorm(n)
    
    dat_pair = list(); dat_pair$prob = rep(1/2, n)
    dat_pair$observe = list(Y1 = Y1, Y2 = Y2, A1 = A1, A2 = A2, X1 = X1, X2 = X2)
    d_pair = C_delta*(Delta(X1)*A1 + Delta(X2)*A2)
    dat_pair$nonnull_ind = d_pair != 0
    dat_pair$nonnull_pos = d_pair > 0
    
    dat_unpair = list(); dat_unpair$prob = rep(1/2, 2*n)
    dat_unpair$observe = list(Y = c(Y1, Y2), A = c(A1, A2), X = rbind(X1, X2))
    d_unpair = C_delta*c(Delta(X1), Delta(X2))
    dat_unpair$nonnull_ind = d_unpair != 0
    dat_unpair$nonnull_pos = d_unpair > 0

    dat = list(dat_pair = dat_pair, dat_unpair = dat_unpair)
  }
  return(dat)
}


