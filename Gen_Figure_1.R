library(CADET)
library(MASS)
library(corrplot)
library(tidyverse)
library(latex2exp)
library(ggsci)
library(ggpubr)

plot_output_dir <- "./plot_output/"
set.seed(2022)
p_val_collection <- rep(NA, times = 1000)
p_naive_collection <- rep(NA, times = 1000)
feat_collection <- rep(NA, times = 1000)
n <- 150
true_clusters <- c(rep(1, 50), rep(2, 50), rep(3, 50))
delta <- 2
q <- 10
mu <- rbind(
  c(delta / 2, rep(0, q - 1)),
  c(rep(0, q - 1), delta / 2), c(rep(0, q - 1), delta / 2)
)
sig <- 1

cov_mat_sim <- matrix(0.4, nrow = q, ncol = q)
diag(cov_mat_sim) <- 1
cov_mat_sim <- cov_mat_sim * sig


png(paste0(plot_output_dir, "covmat_exchange.png"),
  width = 6, height = 6, res = 300, units = "in"
)
corrplot(as.matrix(cov_mat_sim),
  method = "color", type = "full", order = "original",
  is.corr = TRUE,
  tl.col = "black", tl.srt = 45,
  diag = TRUE
)
dev.off()

set.seed(2022)
X <- mvrnorm(n, mu = rep(0, q), Sigma = cov_mat_sim) + mu[true_clusters, ]
X_cluster_concept <- kmeans_estimation(X, k = 2, iter.max = 10, seed = 2021)
final_cluster <- X_cluster_concept$final_cluster
df_feat_1 <- data.frame(final_cluster,
  feat_1 = X[, 1],
  feat_2 = X[, 2]
)
hist_grey <- ggplot(df_feat_1, aes(x = feat_2)) +
  geom_histogram(
    alpha = 0.5,
    aes(y = after_stat(density)),
    position = "identity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  xlab("Feature with no true difference in means") +
  ylab("Density")


hist_color <- ggplot(
  df_feat_1,
  aes(x = feat_2, fill = as.factor(final_cluster))
) +
  geom_histogram(
    alpha = 0.5,
    aes(y = after_stat(density)),
    position = "identity"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  xlab("Feature with no true difference in means") +
  ylab("Density") +
  labs(fill = "") +
  scale_fill_jco(labels = unname(TeX(c("Cluster 1", "Cluster 2"))))



# get a histogram to understand the differences




for (i in c(1:1000)) {
  cat("i=", i, "\n")
  set.seed(i)
  # exchangeable
  X <- mvrnorm(n, mu = rep(0, q), Sigma = cov_mat_sim) + mu[true_clusters, ]
  cluster_num <- sample(c(1:2), 2)
  feat <- sample(c(1:q), 1)
  cl_1_2_inference_demo <- kmeans_inference_1f(X,
    k = 2, cluster_num[1],
    cluster_num[2],
    feat = feat, iso = FALSE,
    sig = NULL, covMat = cov_mat_sim,
    iter.max = 10
  )
  p_val_collection[i] <- cl_1_2_inference_demo$pval
  p_naive_collection[i] <- cl_1_2_inference_demo$p_naive
  feat_collection[i] <- feat
}

df_p0_plot <- data.frame(
  p_selective = (p_val_collection),
  p_naive = (p_naive_collection),
  feature_num = feat_collection
) %>%
  filter(feature_num != 1 & feature_num != 10) %>%
  dplyr::select(-feature_num) %>%
  pivot_longer(everything(), names_to = "p_type", values_to = "p_values") %>%
  mutate(p_type = as.factor(p_type)) %>%
  group_by(p_type) %>%
  mutate(theoretical = ecdf(p_values)(p_values))

# plot ggplot2
type_1_ind <- ggplot(data = df_p0_plot) +
  geom_point(aes(x = theoretical, y = p_values, colour = p_type), size = 0.6) +
  ylab("P-value Quantiles") +
  xlab("Uniform(0,1) Quantiles") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Test") +
  scale_color_jama(labels = unname(TeX(c("$p_{j,naive}$", "$p_{j,k-means}$"))))



# Merging legends
legend_1 <- get_legend(hist_grey)
legend_2 <- get_legend(hist_color)
legend_3 <- get_legend(type_1_ind)
legends <- ggarrange(legend_1, legend_2, legend_3, nrow = 1)

# Combining plots
rm_legend <- function(p) {
  p + theme(legend.position = "none")
}
plots <- ggarrange(rm_legend(hist_grey),
  rm_legend(hist_color),
  rm_legend(type_1_ind),
  align = "hv",
  common.legend = F,
  labels = c("(a)", "(b)", "(c)"),
  ncol = 3, nrow = 1
)

# plots + merged legends
png(paste0(plot_output_dir, "Figure_1.png"),
  height = 4.5,
  width = 14, unit = "in", res = 300
)
ggarrange(plots, legends, heights = c(0.95, 0.05), nrow = 2)
dev.off()
