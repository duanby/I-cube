source("setup.R")

treatment_type = "sparse_pos_bias"
Cd_seq = seq(0, 5, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "methods_unpair",
                        value = c("Crossfit-I-cube", "MaY-I-cube", "linear-BH",
                                  "Crossfit-I-cube-CATE", "MaY-I-cube-RF")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,".Rdata",sep = ""))


treatment_type = "sparse_pos_bias"
Cd_seq = seq(0, 2, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "paired", value = TRUE),
                   list(name = "methods_pair",
                        value = c("pair-Crossfit", "pair-MaY", "unpair-Crossfit", "unpair-MaY")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_pair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/paired_", treatment_type,".Rdata",sep = ""))


treatment_type = "sparse_pos_bias"
eps_seq = seq(0, 0.5, length.out = 6)
result = list()
for (eps in eps_seq) { 
  print(eps)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "paired", value = TRUE),
                   list(name = "eps", value = eps),
                   list(name = "methods_pair",
                        value = c("pair-Crossfit", "pair-MaY", "unpair-Crossfit", "unpair-MaY")),
                   list(name = "C_delta", value = 2))
  result[[as.character(eps)]] = experiment_pair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/paired_", treatment_type,"_mismatch.Rdata",sep = ""))



treatment_type = "sparse_pos_bias"
eps_seq = seq(0, 2.5, length.out = 6)
result = list()
for (eps in eps_seq) { 
  print(eps)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "paired", value = TRUE),
                   list(name = "eps", value = eps),
                   list(name = "methods_pair",
                        value = c("pair-Crossfit", "pair-MaY", "unpair-Crossfit", "unpair-MaY")),
                   list(name = "C_delta", value = 2))
  result[[as.character(eps)]] = experiment_pair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/paired_", treatment_type,"_mismatch_large.Rdata",sep = ""))



treatment_type = "linear"
Cd_seq = seq(0, 5, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "methods_unpair",
                        value = c("Crossfit-I-cube", "MaY-I-cube", "linear-BH")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,".Rdata",sep = ""))


treatment_type = "sparse_oneside"
Cd_seq = seq(0, 5, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "methods_unpair",
                        value = c("Crossfit-I-cube", "MaY-I-cube", "linear-BH")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,".Rdata",sep = ""))


treatment_type = "sparse_twoside"
Cd_seq = seq(0, 5, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "methods_unpair",
                        value = c("Crossfit-I-cube", "MaY-I-cube", "linear-BH")),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_unpair(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,".Rdata",sep = ""))



treatment_type = "subgroup_even"
Cd_seq = seq(0, 5, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "n", value = 2000),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_subgroup(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,".Rdata",sep = ""))


treatment_type = "subgroup_even"
Cd_seq = seq(0, 1, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "n", value = 1000),
                   list(name = "paired", value = TRUE),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_subgroup(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,"_paired.Rdata",sep = ""))



treatment_type = "subgroup_smooth"
Cd_seq = seq(0, 1, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "n", value = 1000),
                   list(name = "paired", value = TRUE),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_subgroup(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,"_paired.Rdata",sep = ""))


treatment_type = "subgroup_sparse"
Cd_seq = seq(0, 1.5, length.out = 6)
result = list()
for (Cd in Cd_seq) { 
  print(Cd)
  para_vary = list(list(name = "treatment_type", value = treatment_type),
                   list(name = "n", value = 1000),
                   list(name = "paired", value = TRUE),
                   list(name = "C_delta", value = Cd))
  result[[as.character(Cd)]] = experiment_subgroup(para_vary)
}
save(result, file=paste(dirname(getwd()),"/result/", treatment_type,"_paired.Rdata",sep = ""))







