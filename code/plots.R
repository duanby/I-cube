source("setup.R")
library(cowplot)
library(ggpubr)
rowSes = function(m) {apply(m, 1, function(x) sd(x, na.rm =  TRUE))/sqrt(rowSums(m, na.rm = TRUE))} #standard error

##################################################################
##################### paired #######################
plot_func = function(result_mat, color_name, exclude_ind = c(), names = NA, se_result_mat = NA,
                     x_label = "scale of treatment effect", y_label,
                     line_type = c("solid", "twodash", "solid", "twodash"),
                     save = FALSE){
  if(mode %in% c("power", "power_pos") & is.na(result_mat[1,1])) {
    result_mat[,1] = 0
    if(!is.na(se_result_mat)) {
      se_result_mat[,1] = 0
    }
  } 
  mode = deparse(substitute(result_mat))
  exclude_methods = rownames(result_mat)[exclude_ind]
  df_power = data.frame(mu_seq = rep(colnames(result_mat), each = nrow(result_mat)),
                        power = as.vector(result_mat),
                        sd = as.vector(se_result_mat),
                        grp = rep(rownames(result_mat),
                                  ncol(result_mat)))
  if (is.na(names)) {
    if (length(exclude_ind)) {
      names = rownames(result_mat)[-exclude_ind]
    } else {
      names = rownames(result_mat)
    }
  }
  p = ggplot(data = subset(df_power, !(grp %in% exclude_methods)),
             aes(x = mu_seq, y = power, group = grp, fill = grp)) +
    geom_line(aes(linetype = grp, color = grp), size = 0.8) +
    geom_point(aes(shape = grp, color = grp), size = 2.5) +
    scale_shape_manual(values = c(21, 22, 21, 22, 24), labels = names) +
    scale_color_manual(values = color_name, labels = names) +
    scale_fill_manual(values = color_name, labels = names) +
    scale_linetype_manual(values = line_type, labels = names) +
    theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "grey"),
        text = element_text(size = 15),
        legend.position = c("none"), legend.text = element_text(size = 10)) +
    xlab(x_label) + ylab(y_label) +
    scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  if (!is.na(se_result_mat)) {
    p = p + geom_errorbar(aes(ymin=power-2*sd, ymax=power+2*sd), width=.05, 
                          position=position_dodge(0.01))
  }
  if (y_label %in% c("FDR", "FDR_neg")) {
    p = p + geom_hline(yintercept = alpha) 
  }
  plot(p)
  return(p)
}


# Figure 3
model = "sparse_pos_bias"
load(file = paste(dirname(getwd()),"/result/", model,".Rdata", sep = ""))
FDR = sapply(result, function(y) {rowMeans(sapply(y, function(x) x$error), na.rm = TRUE)})
FDR_neg = sapply(result, function(y) {rowMeans(sapply(y, function(x) x$error_neg), na.rm = TRUE)})
power = sapply(result, function(y) {rowMeans(sapply(y, function(x) x$power), na.rm = TRUE)})
power_pos = sapply(result, function(y) {rowMeans(sapply(y, function(x) x$power_pos), na.rm = TRUE)})

mode = "FDR"
p = plot_func(result_mat = get(mode), y_label = mode, color_name = c("green3", "red1"),
              exclude_ind = c(2,3,4), names = NA, se_result_mat = NA)
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)

mode = "power_pos"
p = plot_func(result_mat = get(mode), y_label = mode, color_name = c("green3", "red1"),
              exclude_ind = c(2,3,4), names = NA, se_result_mat = NA)
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)



# Figure 4
mode = "FDR"
p = plot_func(result_mat = get(mode), y_label = mode, color_name = c("green3", "deepskyblue2"),
              exclude_ind = c(2,4,5), names = NA, se_result_mat = NA)
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_compare_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)

mode = "FDR_neg"
p = plot_func(result_mat = get(mode), y_label = mode, color_name = c("green3", "deepskyblue2"),
              exclude_ind = c(2,4,5), names = NA, se_result_mat = NA)
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_compare_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)

mode = "power_pos"
p = plot_func(result_mat = get(mode), y_label = mode, color_name = c("green3", "deepskyblue2"),
              exclude_ind = c(2,4,5), names = NA, se_result_mat = NA)
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_compare_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)


# Figure 5
model = "paired_sparse_pos_bias"
load(file = paste(dirname(getwd()),"/result/", model,".Rdata", sep = ""))
power_pos = sapply(result, function(y) {rowMeans(sapply(y, function(x) x$power_pos), na.rm = TRUE)})
mode = "power_pos"
p = plot_func(result_mat = get(mode), y_label = mode,
              color_name = c("tan1","purple","green3", "deepskyblue2"),
              exclude_ind = c(), names = NA, se_result_mat = NA,
              line_type = c("solid","solid", "twodash", "twodash"))
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)

model = "paired_sparse_pos_bias_mismatch"
load(file = paste(dirname(getwd()),"/result/", model,".Rdata", sep = ""))
power_pos = sapply(result, function(y) {rowMeans(sapply(y, function(x) x$power_pos), na.rm = TRUE)})
mode = "power_pos"
p = plot_func(result_mat = get(mode), y_label = mode,
              color_name = c("tan1","purple","green3", "deepskyblue2"),
              exclude_ind = c(), names = NA, se_result_mat = NA,
              line_type = c("solid","solid", "twodash", "twodash"),
              x_label = expression("digree of mismatch"~epsilon))
ggsave(filename = paste(dirname(getwd()),"/figure/",model, "_", mode, ".pdf", sep = ""),
       plot = p, device = "pdf", width = 4, height = 3.6)





