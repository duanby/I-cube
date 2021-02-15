# i-FDR-pair
# dat:       list (Y,A,X) for n subjects
# alpha:     target FDR level
i_FDR_pair = function(dat, alpha, S_model = "RF", alg_type = "Crossfit", cycle_iter = 100){
  n = length(dat$A1)
  random_ind = sample(1:n, 0.5*n) # to not flip the sign
  
  rej_ind1 = rep(TRUE, n)
  rej_ind1[random_ind] = FALSE
  rej_ind1 = eliminate_pair(ini_rej = rej_ind1, alpha = alpha/2, dat = dat, 
                       S_model = S_model, alg_type = alg_type, cycle_iter = cycle_iter)
  rej_ind2 = rep(TRUE, n)
  rej_ind2[-random_ind] = FALSE
  rej_ind2 = eliminate_pair(ini_rej = rej_ind2, alpha = alpha/2, dat = dat,
                       S_model = S_model, alg_type = alg_type, cycle_iter = cycle_iter)
  final_rej = (rej_ind1 + rej_ind2 > 0)
  return(final_rej)
}


eliminate_pair = function(ini_rej, alpha, dat, S_model, alg_type, cycle_iter){
  n = length(dat$A1); df = data.frame(dat)
  
  rej_ind = ini_rej; iter = 0; fdr_est = 1
  incre = (dat$A1 - dat$A2)*(dat$Y1 - dat$Y2); sign = factor(incre <= 0)
  
  if(S_model == "RF" & alg_type == "Crossfit") {
    pred_dat = data.frame(sign = sign, dat$X1, dat$X2, dat$Y1, dat$Y2)
  }
  
  while (fdr_est > alpha & any(rej_ind)) {
    if (iter %% cycle_iter == 0) {
        if (alg_type == "Crossfit") {
          if (S_model == "CATE") {
            est_effect = incre
            effect_rf = randomForest(y = est_effect[!rej_ind],
                                     x = cbind(dat$X1[!rej_ind,], dat$X2[!rej_ind,],
                                               dat$Y1[!rej_ind], dat$Y2[!rej_ind]))
            S = predict(effect_rf, newdata = cbind(dat$X1, dat$X2, dat$Y1, dat$Y2)); S[!rej_ind] = NA
          } else if (S_model == "RF") {
            sign_model = randomForest(sign ~ ., data = pred_dat[!rej_ind,])
            S = predict(sign_model, newdata = pred_dat, type = "prob")[,1]; S[!rej_ind] = NA
          }
        } else if (alg_type == "MaY") {
          if (S_model == "CATE") {
            est_effect = incre
            effect_rf = randomForest(y = est_effect[!rej_ind],
                                     x = cbind(dat$X1[!rej_ind,], dat$X2[!rej_ind,]))
            S = predict(effect_rf, newdata = cbind(dat$X1, dat$X2)); S[!rej_ind] = NA
          } else if (S_model %in% c("RF")) {
            rf_y1 = randomForest(Y1 ~ . - A1 - A2 - Y2, data = df[!rej_ind,]); y1_hat = predict(rf_y1, newdata = df)
            rf_y2 = randomForest(Y2 ~ . - A1 - A2 - Y1, data = df[!rej_ind,]); y2_hat = predict(rf_y2, newdata = df)
            pred_dat = data.frame(sign = sign, dat$X1, dat$X2,
                                  y1 = y1_hat*rej_ind + dat$Y1*(!rej_ind),
                                  y2 = y2_hat*rej_ind + dat$Y2*(!rej_ind))
            sign_model = randomForest(sign ~ ., data = pred_dat[!rej_ind,])
            S = predict(sign_model, newdata = pred_dat, type = "prob")[,1]; S[!rej_ind] = NA
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



