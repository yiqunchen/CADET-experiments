library(CADET)
library(MASS)
library(ggplot2)
library(latex2exp)
library(ggpubr)

plot_output_dir <- "./plot_output/"
set.seed(2022)
n <- 30
true_clusters <- c(rep(1, 10), rep(2, 10), rep(3, 10))
delta <- 5
q <- 2
mu <- rbind(
  c(0, rep(delta / 2, q - 1)),
  c(rep(delta / 2, q - 1), 0),
  c(-delta / 2, -rep(delta / 2, q - 1))
)
sig <- 0.2

cov_mat_sim <- matrix(0.4, nrow = q, ncol = q)
diag(cov_mat_sim) <- 1
cov_mat_sim <- cov_mat_sim * sig

X <- mvrnorm(n, mu = rep(0, q), Sigma = cov_mat_sim) + mu[true_clusters, ]
cluster_num <- c(1, 2)
cluster_1 <- cluster_num[1]
cluster_2 <- cluster_num[2]
feat <- 1
real_diff_in_means <- mean(mu[cluster_num[1], feat]) - mean(mu[cluster_num[2], feat])

estimated_k_means <- kmeans_estimation(X, k = 3, iter.max = 10, seed = 2022)
estimated_final_cluster <- estimated_k_means$cluster[[estimated_k_means$iter]]
all_T_clusters <- do.call(rbind, estimated_k_means$cluster)
all_T_centroids <- estimated_k_means$centers
T_length <- nrow(all_T_clusters)

cl_1_2_inference <- kmeans_inference_1f(X,
  k = 3, cluster_num[1],
  cluster_num[2],
  feat = feat,
  iso = FALSE,
  sig = NULL,
  seed = 2022,
  covMat = cov_mat_sim
)

# construct contrast vector
v_vec <- rep(0, times = nrow(X))
mean_indicator_vec <- (estimated_final_cluster == cluster_1) | (estimated_final_cluster == cluster_2)
a_vec[mean_indicator_vec] <- 1 / (sum(mean_indicator_vec))
v_vec[estimated_final_cluster == cluster_1] <- 1 / (sum(estimated_final_cluster == cluster_1))
v_vec[estimated_final_cluster == cluster_2] <- -1 / (sum(estimated_final_cluster == cluster_2))

n1 <- sum(estimated_final_cluster == cluster_1)
n2 <- sum(estimated_final_cluster == cluster_2)
squared_norm_nu <- 1 / n1 + 1 / n2
v_norm <- sqrt(squared_norm_nu) # recycle this computed value
# compute XTv
diff_means_feat <- mean(X[estimated_final_cluster == cluster_1, feat]) -
  mean(X[estimated_final_cluster == cluster_2, feat])

colMeans(X[estimated_final_cluster == cluster_1, ])
colMeans(X[estimated_final_cluster == cluster_2, ])

p_orig <- ggplot(data.frame(X), aes(
  x = X1, y = X2,
  col = as.factor(estimated_final_cluster)
)) +
  geom_point(cex = 2) +
  xlab("Feature 1") +
  ylab("Feature 2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  xlim(c(-4, 5)) +
  ylim(c(-6, 4)) +
  scale_colour_manual(values = c("dodgerblue3", "rosybrown", "orange")) +
  ggtitle(TeX(
    paste0("Original data ($\\phi$=$X_j^T \\nu$=$", round(diff_means_feat, 2), ")")
  )) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_vline(
    xintercept = mean(X[estimated_final_cluster == cluster_1, 1]),
    col = "dodgerblue3", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(X[estimated_final_cluster == cluster_2, 1]),
    col = "rosybrown", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = mean(X[estimated_final_cluster == cluster_1, 2]),
    col = "dodgerblue3", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = mean(X[estimated_final_cluster == cluster_2, 2]),
    col = "rosybrown", linetype = "dashed"
  )


# ok happy with the initial plot
# perturbation func
get_perturbed_x <- function(x, phi, v, diff_in_means, feat, Sigma) {
  x_phi <- x + as.matrix((phi - diff_in_means) * v / (sum(v^2))) %*% t(Sigma[feat, ]) / Sigma[feat, feat]
  return(x_phi)
}

test_recover <- get_perturbed_x(X,
  phi = diff_means_feat, v = v_vec,
  diff_in_means = diff_means_feat,
  feat = feat, Sigma = cov_mat_sim
)
# this should be true
all.equal(test_recover, X)



X_move_larger <- get_perturbed_x(X,
  phi = 0, v = v_vec,
  diff_in_means = diff_means_feat,
  feat = feat, Sigma = cov_mat_sim
)

p_large <- ggplot(data.frame(X_move_larger), aes(
  x = X1, y = X2,
  col = as.factor(estimated_final_cluster)
)) +
  geom_point(cex = 2) +
  xlab("Feature 1") +
  ylab("Feature 2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  xlim(c(-4, 5)) +
  ylim(c(-6, 4)) +
  scale_colour_manual(values = c("dodgerblue3", "rosybrown", "orange")) +
  ggtitle(TeX(
    paste0("Perturbed data ($\\phi$=", round((0), 2), ")")
  )) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_vline(
    xintercept = mean(X_move_larger[estimated_final_cluster == cluster_1, 1]),
    col = "dodgerblue3", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(X_move_larger[estimated_final_cluster == cluster_2, 1]),
    col = "rosybrown", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = mean(X_move_larger[estimated_final_cluster == cluster_1, 2]),
    col = "dodgerblue3", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = mean(X_move_larger[estimated_final_cluster == cluster_2, 2]),
    col = "rosybrown", linetype = "dashed"
  )


X_move_smaller <- get_perturbed_x(X,
  phi = -6, v = v_vec,
  diff_in_means = diff_means_feat,
  feat = feat, Sigma = cov_mat_sim
)

p_small <- ggplot(data.frame(X_move_smaller), aes(
  x = X1, y = X2,
  col = as.factor(estimated_final_cluster)
)) +
  geom_point(cex = 2) +
  xlab("Feature 1") +
  ylab("Feature 2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none") +
  xlim(c(-4, 5)) +
  ylim(c(-6, 4)) +
  scale_colour_manual(values = c("dodgerblue3", "rosybrown", "orange")) +
  ggtitle(TeX(
    paste0("Perturbed data ($\\phi$=", round((-5)), ")")
  )) +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_vline(
    xintercept = mean(X_move_smaller[estimated_final_cluster == cluster_1, 1]),
    col = "dodgerblue3", linetype = "dashed"
  ) +
  geom_vline(
    xintercept = mean(X_move_smaller[estimated_final_cluster == cluster_2, 1]),
    col = "rosybrown", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = mean(X_move_smaller[estimated_final_cluster == cluster_1, 2]),
    col = "dodgerblue3", linetype = "dashed"
  ) +
  geom_hline(
    yintercept = mean(X_move_smaller[estimated_final_cluster == cluster_2, 2]),
    col = "rosybrown", linetype = "dashed"
  )


phi_seq <- seq(-5, 5, by = 0.1)
delta_matrix <- matrix(nrow = length(phi_seq), ncol = 2)
for (i in seq_along(phi_seq)) {
  x_1 <- get_perturbed_x(X,
    phi = phi_seq[i], v = v_vec,
    diff_in_means = diff_means_feat,
    feat = feat, Sigma = cov_mat_sim
  )
  delta_matrix[i, ] <- colMeans(x_1[estimated_final_cluster == cluster_1, ]) -
    colMeans(x_1[estimated_final_cluster == cluster_2, ])
}

panel_d <- data.frame(
  phi = phi_seq,
  delta_1 = delta_matrix[, 1],
  delta_2 = delta_matrix[, 2]
) %>%
  pivot_longer(cols = c(delta_1, delta_2))


p2_d <- panel_d %>% ggplot(aes(x = phi, y = value, color = name)) +
  geom_line(linewidth = 1.5) +
  xlab(TeX("$\\phi$")) +
  ylab("Difference in feature means") +
  theme_classic(base_size = 14) +
  scale_color_discrete("", labels = c("Feature 1", "Feature 2")) +
  theme(legend.position = c(0.8, 0.2))


plots_fig2 <- ggarrange(p_orig,
  p_large,
  p_small,
  p2_d,
  align = "hv",
  common.legend = F,
  labels = c("(a)", "(b)", "(c)", "(d)"),
  ncol = 4, nrow = 1
)

png(paste0(plot_output_dir, "Figure_2.png"),
  height = 4,
  width = 14.5, unit = "in", res = 300
)
plots_fig2
dev.off()
