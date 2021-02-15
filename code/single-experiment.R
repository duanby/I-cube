experiment_unpair = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    rejections = list()
    
    dat = sample_generator(n = n, C_delta = C_delta, 
                           treatment_type = treatment_type, control_type = control_type,
                           m = m, level = level, paired = FALSE, eps = NA)

    if ("I-cube-bet" %in% methods_unpair) {
      rejections[["I-cube-bet"]] = 
        i_FDR_unpair(dat = dat$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "no_split", S_model = "linear")
    }
    if ("Crossfit-I-cube" %in% methods_unpair) {
      rejections[["Crossfit-I-cube"]] = 
        i_FDR_unpair(dat = dat$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "Crossfit", S_model = "RF")
    }
    if ("Crossfit-I-cube-CATE" %in% methods_unpair) {
      rejections[["Crossfit-I-cube-CATE"]] = 
        i_FDR_unpair(dat = dat$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "Crossfit", S_model = "CATE")
    }
    if ("MaY-I-cube" %in% methods_unpair) {
      rejections[["MaY-I-cube"]] = 
        i_FDR_unpair(dat = dat$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "MaY", S_model = "CATE")
    }
    if ("MaY-I-cube-RF" %in% methods_unpair) {
      rejections[["MaY-I-cube-RF"]] = 
        i_FDR_unpair(dat = dat$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "MaY", S_model = "RF")
    }
    if ("linear-BH" %in% methods_unpair) {
      rejections[["linear-BH"]] = 
        linear_FDR(dat = dat$observe, alpha = alpha)
    }
    
    power = sapply(rejections, function(x) {sum(dat$nonnull_ind & x)/sum(dat$nonnull_ind)})
    error = sapply(rejections, function(x) {sum((!dat$nonnull_ind) & x)/max(sum(x), 1)})
    
    power_pos = sapply(rejections, function(x) {sum(dat$nonnull_pos & x)/sum(dat$nonnull_pos)})
    error_neg = sapply(rejections, function(x) { sum((!dat$nonnull_pos) & x)/max(sum(x), 1) })

    return(list(power = power, error = error,
                error_neg = error_neg, power_pos = power_pos))
  }
  rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  # rejections = lapply(1:R, wrapper_func)
  return(rejections)
}


experiment_pair = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)
    
    dat = sample_generator(n = n, C_delta = C_delta, 
                           treatment_type = treatment_type, control_type = control_type,
                           m = m, level = level, paired = TRUE,  eps = eps)
    rejections_pair = list()
    if ("pair-Crossfit" %in% methods_pair) {
      rejections_pair[["pair-Crossfit"]] = 
        i_FDR_pair(dat = dat$dat_pair$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "Crossfit", S_model = "RF")
    }
    if ("pair-MaY" %in% methods_pair) {
      rejections_pair[["pair-MaY"]] = 
        i_FDR_pair(dat = dat$dat_pair$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "MaY", S_model = "CATE")
    }
    power_pair = sapply(rejections_pair, function(x)
      {sum(dat$dat_pair$nonnull_ind & x)/sum(dat$dat_pair$nonnull_ind)})
    error_pair = sapply(rejections_pair, function(x)
      {sum((!dat$dat_pair$nonnull_ind) & x)/max(sum(x), 1)})
    
    power_pos_pair = sapply(rejections_pair, function(x) 
      {sum(dat$dat_pair$nonnull_pos & x)/sum(dat$dat_pair$nonnull_pos)})
    error_neg_pair = sapply(rejections_pair, function(x)
      { sum((!dat$dat_pair$nonnull_pos) & x)/max(sum(x), 1) })
    
    rejections_unpair = list()
    if ("unpair-Crossfit" %in% methods_pair) {
      rejections_unpair[["unpair-Crossfit"]] = 
        i_FDR_unpair(dat = dat$dat_unpair$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "Crossfit", S_model = "RF")
    }
    if ("unpair-MaY" %in% methods_pair) {
      rejections_unpair[["unpair-MaY"]] = 
        i_FDR_unpair(dat = dat$dat_unpair$observe, alpha = alpha, cycle_iter = cycle_iter,
                     alg_type = "MaY", S_model = "CATE")
    }
    
    power_unpair = sapply(rejections_unpair, function(x)
      {sum(dat$dat_unpair$nonnull_ind & x)/sum(dat$dat_unpair$nonnull_ind)})
    error_unpair = sapply(rejections_unpair, function(x)
      {sum((!dat$dat_unpair$nonnull_ind) & x)/max(sum(x), 1)})
    
    power_pos_unpair = sapply(rejections_unpair, function(x)
      {sum(dat$dat_unpair$nonnull_pos & x)/sum(dat$dat_unpair$nonnull_pos)})
    error_neg_unpair = sapply(rejections_unpair, function(x)
      {sum((!dat$dat_unpair$nonnull_pos) & x)/max(sum(x), 1) })
    
    return(list(power = c(power_pair, power_unpair), error = c(error_pair, error_unpair),
                error_neg = c(error_neg_pair, error_neg_unpair),
                power_pos = c(power_pos_pair, power_pos_unpair)))
  }
  rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  # rejections = lapply(1:R, wrapper_func)
  return(rejections)
}

experiment_subgroup = function(para_vary){
  source("input.R", local = TRUE)
  for (single_para_vary in para_vary) {
    assign(single_para_vary$name, single_para_vary$value)
  }
  
  wrapper_func = function(i){
    print(i)

    dat = sample_generator(n = n, C_delta = C_delta, 
                           treatment_type = treatment_type, control_type = control_type,
                           m = m, level = level, paired = paired, eps = NA)
    if (paired) {dat = dat$dat_pair}
    result = 
      subgroup_FDR(dat$observe, alpha = alpha, nonnull_ind = dat$nonnull_ind, paired = paired)
    
    return(result)
  }
  rejections = mclapply(1:R, wrapper_func, mc.cores = detectCores())
  # rejections = lapply(1:R, wrapper_func)
  return(rejections)
}




