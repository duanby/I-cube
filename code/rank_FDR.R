# i-FDR-unpair
# dat:       list (Y,A,X) for n subjects
# alpha:     target FDR level
i_FDR_unpair = function(dat, alpha, alg_type = "Crossfit", S_model = "RF", cycle_iter = 100){
  n = length(dat$A)
  
  if (alg_type == "no_split") {
    rej_ind = rep(TRUE, n)
    # random_ind = sample(1:n, 0.1*n); rej_ind[random_ind] = FALSE
    final_rej = eliminate(ini_rej = rej_ind, alpha = alpha, dat = dat,
                         S_model = S_model, alg_type = alg_type, cycle_iter = cycle_iter)
  } else {
    random_ind = sample(1:n, 0.5*n) 
    
    rej_ind1 = rep(TRUE, n)
    rej_ind1[random_ind] = FALSE
    rej_ind1 = eliminate(ini_rej = rej_ind1, alpha = alpha/2, dat = dat,
                         S_model = S_model, alg_type = alg_type, cycle_iter = cycle_iter)
    rej_ind2 = rep(TRUE, n)
    rej_ind2[-random_ind] = FALSE
    rej_ind2 = eliminate(ini_rej = rej_ind2, alpha = alpha/2, dat = dat, 
                         S_model = S_model, alg_type = alg_type, cycle_iter = cycle_iter)
    final_rej = (rej_ind1 + rej_ind2 > 0)
  }
  return(final_rej)
}


eliminate = function(ini_rej, alpha, dat, S_model, alg_type, cycle_iter){
  n = length(dat$A); df = data.frame(dat)
  
  prob_bd = 1/2
  rej_ind = ini_rej; iter = 0; fdr_est = 1
  #prepare: increment
  if (alg_type == "no_split"){
    rf_y = lm(Y ~ (. - A)^2, data = df)
    # rf_y = randomForest(Y ~ . - A, data = df)
    yhat = predict(rf_y)
    e = dat$Y - yhat
    incre = 2*(2*dat$A - 1)*e; sign = factor(incre <= 0)
  } else if (alg_type == "Crossfit"){
    if (length(unique(dat$Y)) == 2) {
      df$Y = as.factor(dat$Y)
      rf_y = randomForest(Y ~ . - A, data = df)
      yhat = predict(rf_y, type = "prob")[,2]
    } else {
      rf_y = randomForest(Y ~ . - A, data = df)
      yhat = predict(rf_y)
    }
    e = dat$Y - yhat
    incre = 2*(2*dat$A - 1)*e; sign = factor(incre <= 0)
    if(S_model == "RF") {pred_dat = data.frame(sign = sign, e = e, y = dat$Y, dat$X)}
  } else {
    if(length(unique(dat$Y)) == 2) {
      df$Y = as.factor(dat$Y)
      rf_y = randomForest(Y ~ . - A, data = df[!rej_ind,])
      yhat = predict(rf_y, newdata = df, type = "prob")[,2]
    } else {
      rf_y = randomForest(Y ~ . - A, data = df[!rej_ind,])
      yhat = predict(rf_y, newdata = df)
    }
    e = dat$Y - yhat; 
    incre = 2*(2*dat$A - 1)*e; sign = factor(incre <= 0)
    df_e = data.frame(e = e, ind_e = e > 0, dat$X)
  }
  
  while (fdr_est > alpha & any(rej_ind)) {
    if (iter %% cycle_iter == 0) { ## update S score every 100 iterations
      if (alg_type == "no_split") {
        S = em_bet(y = incre, x = dat$X, reveal_ind = !rej_ind)
      } else if (alg_type == "Crossfit") {
        if (S_model == "RF") {
          sign_model = randomForest(sign ~ ., data = pred_dat[!rej_ind,])
          S = predict(sign_model, newdata = pred_dat, type = "prob")[,1]; S[!rej_ind] = NA
        } else if (S_model == "CATE") {
          potential_rf = randomForest(Y ~ ., data = df[!rej_ind,])
          est_true = predict(potential_rf)
          pseudo_treat = df[!rej_ind,] %>% mutate(A = rep(1, sum(!rej_ind)))
          pseudo_control = df[!rej_ind,] %>% mutate(A = rep(0, sum(!rej_ind)))
          est_treat = predict(potential_rf, newdata = pseudo_treat); est_control = predict(potential_rf, newdata = pseudo_control)
          
          est_effect = 4*(dat$A[!rej_ind] - 1/2)*
            (dat$Y[!rej_ind] - est_true) + est_treat - est_control
          est_effect_full = rep(NA, n); est_effect_full[!rej_ind] = est_effect
          df_est = data.frame(y = est_effect_full, x = dat$X)
          effect_rf = randomForest(y ~., data = df_est[!rej_ind,])
          S = predict(effect_rf, newdata = df_est); S[!rej_ind] = NA
        }
      } else if (alg_type == "MaY") {
        if (S_model == "RF") {
          rf_e = randomForest(e ~ ., data = df_e[!rej_ind,]); ehat = predict(rf_e, newdata = df_e)
          pred_dat = data.frame(sign = sign, e = ehat, y = yhat*rej_ind + dat$Y*(!rej_ind), dat$X)
          sign_model = randomForest(sign ~ ., data = pred_dat[!rej_ind,])
          S = predict(sign_model, newdata = pred_dat, type = "prob")[,1]; S[!rej_ind] = NA
        } else if (S_model == "CATE") {
          potential_rf = randomForest(Y ~ ., data = df[!rej_ind,])
          est_true = predict(potential_rf)
          pseudo_treat = df[!rej_ind,] %>% mutate(A = rep(1, sum(!rej_ind)))
          pseudo_control = df[!rej_ind,] %>% mutate(A = rep(0, sum(!rej_ind)))
          est_treat = predict(potential_rf, newdata = pseudo_treat); est_control = predict(potential_rf, newdata = pseudo_control)
          
          est_effect = 4*(dat$A[!rej_ind] - 1/2)*
            (dat$Y[!rej_ind] - est_true) + est_treat - est_control
          effect_rf = randomForest(y = est_effect, x = dat$X[!rej_ind,])
          S = predict(effect_rf, newdata = dat$X); S[!rej_ind] = NA
        }
      }
    }
    rej_ind[which.min(S)] = FALSE
    S[which.min(S)] = NA
    
    r_neg = sum(incre[rej_ind] <= 0); r_pos = sum(incre[rej_ind] > 0)
    fdr_est = (r_neg + 1)/max(r_pos,1)

    iter = iter + 1
  }
  
  final_rej = rej_ind; final_rej[incre <= 0] = FALSE
  return(final_rej)
}

em_bet <- function(y, x, reveal_ind, iter = 10){
  y_obs = abs(y)*(!reveal_ind) + y*reveal_ind
  rank_deficient = any(is.na(cor(model.matrix(~ .^2 - 1, data.frame(x))[reveal_ind,])))
  if (sum(reveal_ind) == 0 | rank_deficient) {
    mu = rep(0.1, length(y_obs))
  } else {
    ini_df = data.frame(Y = y_obs, x = x) #close to corret modeling
    ini_regress = lm(Y ~ .^2, data = ini_df[reveal_ind,])
    mu = predict(ini_regress, newdata = ini_df)
  }
  for (i in 1:iter) {
    w = 1/(1 + exp(-2*mu*y_obs))
    w[reveal_ind] <- 1 
    weight_df = data.frame(Y = (2*w - 1)*y_obs, x = x)
    mu = lm(Y ~ .^2, data = weight_df)$fitted.values #not correctly learned
  }
  w[reveal_ind] = NA
  
  # if (sum(reveal_ind) == 0 & sum(w < 0.1) > sum(w > 0.9)) { #bet on positive
  #   w = 1 - w
  # }
  return(w)
}


#############################################
############### linear-BH ###################
#############################################
linear_FDR = function(dat, alpha) {
  df = data.frame(Y = dat$Y, X = dat$X)
  Y_model_treat = lm(Y ~ .^2, data = df[dat$A == 1,])
  Y_model_control = lm(Y ~ .^2, data = df[dat$A == 0,])
  
  pred_treat = predict(Y_model_treat, newdata = df, se.fit = TRUE)
  pred_control = predict(Y_model_control, newdata = df, se.fit = TRUE)
  Y_treat = dat$Y*(dat$A == 1) + pred_treat$fit*(dat$A == 0)
  Y_control = pred_control$fit*(dat$A == 1) + dat$Y*(dat$A == 0)
  var = (pred_control$se.fit^2 + pred_treat$residual.scale^2)*(dat$A == 1) +
    (pred_treat$se.fit^2 + pred_control$residual.scale^2)*(dat$A == 0)
  p_val = 1 - pnorm((Y_treat - Y_control)/sqrt(var))
  
  rejection = p.adjust(p = p_val, method = "BH") < alpha
  return(rejection)
}



#############################################
############### subgroup ###################
#############################################
subgroup_FDR = function(dat, alpha, nonnull_ind, paired){
  rejections = list()
  df = data.frame(dat)
  if(!paired) {
    grp <- predict(partykit::ctree(Y ~ . - A, data = df))
    label_grp = unique(grp); n_grp = length(label_grp)
    e = dat$Y - randomForest(Y ~ . - A, data = df)$predicted
    p = rep(0, length(n_grp)); x = matrix(nrow = n_grp, ncol = ncol(dat$X) + 1)
    nonnull_grp = rep(FALSE, n_grp)
    for(i in 1:n_grp){
      p[i] <- wilcox.test((2*dat$A[grp == label_grp[i]] - 1)*e[grp == label_grp[i]],
                          rep(0, sum(grp == label_grp[i])),
                          paired = TRUE, alternative = "greater")$p.value 
      x[i,] <- c(colMeans(dat$X[grp == label_grp[i],]), mean(dat$Y[grp == label_grp[i]]))
      nonnull_grp[i] = any(nonnull_ind[grp == label_grp[i]])
    }
  } else {
    diff_y = dat$Y1 - dat$Y2; diff_a = dat$A1 - dat$A2
    grp = dat$X1[,1]*(dat$X1[,2] == 1) + (1000 + dat$X1[,1])*(dat$X1[,2] == 2) 
    label_grp = unique(grp); n_grp = length(label_grp)
    p = rep(0, length(n_grp)); x = matrix(nrow = n_grp, ncol = 2*ncol(dat$X1) + 2)
    nonnull_grp = rep(FALSE, n_grp)
    for(i in 1:n_grp){
      p[i] <- wilcox.test(diff_a[grp == label_grp[i]]*diff_y[grp == label_grp[i]],
                          rep(0, sum(grp == label_grp[i])),
                          paired = TRUE, alternative = "greater")$p.value
      if(sum(grp == label_grp[i]) == 1) {
        x[i,] <- c(dat$X1[grp == label_grp[i],], dat$X2[grp == label_grp[i],],
                   dat$Y1[grp == label_grp[i]], dat$Y2[grp == label_grp[i]])
      } else {
        x[i,] <- c(colMeans(cbind(dat$X1[grp == label_grp[i],], dat$X2[grp == label_grp[i],])),
                   mean(dat$Y1[grp == label_grp[i]]), mean(dat$Y2[grp == label_grp[i]])) 
      }
      nonnull_grp[i] = any(nonnull_ind[grp == label_grp[i]])
    }
  }
  
  rejections[["BH"]] = p.adjust(p = p, method = "BH") < alpha
  
  rejections[["adaptive"]] = i_FDR(P = p, x = x, alpha = alpha, S_model = FALSE,
                                   mask_fun = "tent", mask_para = 0.5, structure = "sequence")
  
  rejections[["interactive"]] = i_FDR(P = p, x = x, alpha = alpha,
                                     mask_fun = "tent", mask_para = 0.5, structure = "sequence")
  
  power = sapply(rejections, function(x) {sum(nonnull_grp & x)/sum(nonnull_grp)})
  error = sapply(rejections, function(x) {sum((!nonnull_grp) & x)/max(sum(x), 1)})
  return(list(power = power, error = error, total_n = n_grp))
}



# i-FDR
# P:         vector of p-values of length n
# x:         side information 
# alpha:     target FWER level
# mask_fun:  name of the masking function
# mask_para: a vector of parameters in the masking function 
#            eg. [p_l, p_u] for gap function
i_FDR = function(P, x, alpha, mask_fun, mask_para, structure, S_model = TRUE,
                 d = 5, delta = 0.05, tree_obj = NULL){
  # generate masked p-values and missing bits
  h_g = mask_gen(P, mask_fun, mask_para); h = h_g$h; g = h_g$g
  if(mask_fun %in% c("gap", "gap-railway")) {mask_ind = !(h == 0)}
  
  
  n = length(P)
  rej_ind = rep(TRUE, n)  ## start with all the hypotheses in the rejection set
  r_neg = sum(h[rej_ind] == -1)
  fdr_est = (r_neg + 1)/max(sum(rej_ind) - r_neg, 1)
  iter = 0
  while (fdr_est > alpha & any(rej_ind)) {
    masked_P = g*rej_ind + P*(1 - rej_ind)
    if (S_model) { 
      if (iter %% 100 == 0) { ## update S score every 100 iterations
        if (mask_fun %in% c("tent", "railway")) {
          S = em_mixture(masked_P, x, rej_ind, mask_fun, mask_para, structure)
        } else {
          S = em_mixture(masked_P, x, rej_ind & mask_ind, mask_fun, mask_para, structure)
        }
      }
    } else {
      S = -masked_P
    }
    
    if(structure == "grid") {
      rej_ind = exclu_search(S, x, rej_ind, d, delta)
    } else if (structure == "tree") {
      rej_ind = exclu_search_tree(S, tree_obj, rej_ind)
    } else if (structure == "sequence") {
      if(sum(rej_ind) > 10) {
        cut_S = quantile(S[rej_ind], probs = c(0.1))
      } else {
        cut_S = min(S[rej_ind])
      }
      rej_ind[S <= cut_S] = FALSE
    }
    
    r_neg = sum(h[rej_ind] == -1)
    fdr_est = (r_neg + 1)/max(sum(rej_ind) - r_neg, 1)
    iter = iter + 1
  }
  rejections = (rej_ind & h == 1)
  return(rejections)
}

# masked p-value generator
# P:         vector of p-values of length n
# mask_fun:  name of the masking function
# mask_para: a vector of parameters in the masking function
mask_gen = function(P, mask_fun, mask_para){
  if (mask_fun == "tent"){
    h <- 1*(P < mask_para) + (-1)*(P >= mask_para)
    g <- pmin(P, mask_para/(1 - mask_para)*(1 - P))
  } else if (mask_fun == "railway"){
    h <- 1*(P < mask_para) + (-1)*(P >= mask_para)
    g <- P*(P < mask_para) +
      mask_para/(1 - mask_para)*(P - mask_para)*(P >= mask_para) 
  } else if (mask_fun == "gap"){
    h <- 1*(P < mask_para[1]) + (-1)*(P > mask_para[2])
    g <- P*(P < mask_para[1]) +
      P*(P >= mask_para[1] & P <= mask_para[2]) +
      mask_para[1]/(1 - mask_para[2])*(1 - P)*(P > mask_para[2])
  }
  else if (mask_fun == "gap-railway"){
    h <- 1*(P < mask_para[1]) + (-1)*(P > mask_para[2])
    g <- P*(P < mask_para[1]) +
      P*(P >= mask_para[1] & P <= mask_para[2]) +
      mask_para[1]/(1 - mask_para[2])*(P - mask_para[2])*(P > mask_para[2]) #gap-railway
  }
  return(list(h = h, g = g))
}



# EM algorithm under mixture model to get the non-null likelihood score
# masked_P:  vector of masked p-value information
# x:         side information 
# rej_ind:   vector of indicators for whether a hypothesis is in the rejection set
# mask_fun:  name of the masking function
# mask_para: a vector of parameters in the masking function
# iter:      number of EM iterations, default to five
em_mixture = function(masked_P, x, rej_ind, mask_fun, mask_para, structure, iter = 5, df = 3){
  masked_P[masked_P < 10^(-10)] = 10^(-10)  ## avoid Inf in calculation
  masked_Z = qnorm(1 - masked_P) ## translate p-values into Z-scores
  if (mask_fun == "tent") {
    inv_Z = qnorm((1 - mask_para)/mask_para*(1 - pnorm(masked_Z)))
  } else if (mask_fun == "railway") {
    inv_Z = qnorm((1 - mask_para)/mask_para*(pnorm(masked_Z) - 1 + mask_para)) 
  } else if (mask_fun == "gap") {
    inv_Z = qnorm((1 - mask_para[2])/mask_para[1]*(1 - pnorm(masked_Z)))
  }else if (mask_fun == "gap-railway") {
    inv_Z = qnorm((1 - mask_para[2])/mask_para[1]*(pnorm(masked_Z) - 1 + mask_para[1]))
  }
  inv_Z[!rej_ind] = 0
  
  pi_set = rep(0.1, length(masked_P)); mu = 1 ## initial values for parameters in the EM algorithm
  for (i in 1:iter) {
    mu = mask_para %>% ifelse(mask_fun %in% c("tent", "railway"), ., .[1]) %>%
      max(mu, qnorm(1 - .) + 0.5)
    # likelihood of each case with respect to the latent labels w and q
    # a = pi_set*dnorm(masked_Z - mu); b = (1 - pi_set)*dnorm(masked_Z)
    # c = pi_set*dnorm(inv_Z - mu); d = (1 - pi_set)*dnorm(inv_Z)
    a = pi_set*dnorm(masked_Z - mu); b = (1 - pi_set)*dnorm(masked_Z)
    if (mask_fun %in% c("gap", "gap-railway")) {
      drv = dnorm(masked_Z)/dnorm(inv_Z)*(1 - mask_para[2])/mask_para[1]
    } else if (mask_fun %in% c("tent", "railway")) {
      drv = dnorm(masked_Z)/dnorm(inv_Z)*(1 - mask_para)/mask_para
    }
    #c = pi_set*dnorm(inv_Z - mu)*drv; d = (1 - pi_set)*dnorm(inv_Z)*drv
    
    a = ifelse(rej_ind, a, a/(a+b)); a[is.na(a)] = 0.5
    b = ifelse(rej_ind, b, 1 - a)
    c = ifelse(rej_ind, pi_set*dnorm(inv_Z - mu)*drv, 0)
    d = ifelse(rej_ind, (1 - pi_set)*dnorm(inv_Z)*drv, 0)
    sum_abcd = a + b + c + d
    q = (a+c)/sum_abcd
    a = a/sum_abcd; b = b/sum_abcd; c = c/sum_abcd; d = d/sum_abcd
    # update latent labels: w, q
    # w = ifelse(rej_ind, 1/(1 + (c + d)/(a + b)), 1)
    # q = ifelse(rej_ind, 1/(1 + (b + d)/(a + c)), 1/(1 + b/a))
    # update parameters: pi_set, mu
    if (structure == "grid") {
      phi_x = bs(x[,1], df = df); phi_y = bs(x[,2], df); 
      phi = phi_x[,rep(1:df, each = df)] * phi_y[,rep(1:df, times = df)]
      pi_set = glm(q~phi, family = quasibinomial())$fitted.values
    } else if (structure == "tree") {
      pi_set = activeSet(isomat = x[,c(2,1)], mySolver = "LS", y = q,
                         weights = rep(1, length(masked_P)))$x
    } else if (structure == "sequence") {
      phi_x = bs(x[,1], df = df); phi_y = bs(x[,2], df); phi_z = bs(x[,3], df = df)
      phi = phi_x*phi_y*phi_z
      pi_set = glm(q~phi, family = quasibinomial())$fitted.values
      # temp = a + c
      # pi_set = glm(temp ~ phi, family = quasibinomial())$fitted.values
    }
    # mu = sum(q*w*masked_Z + q*(1 - w)*inv_Z)/sum(q)
    mu = sum(a*masked_Z + c*inv_Z)/sum(a + c)
  }
  return(q)
}

