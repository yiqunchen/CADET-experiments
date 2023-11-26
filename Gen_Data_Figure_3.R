set.seed(2022)
library(CADET)
library(MASS)
library(corrplot)
library(tidyverse)
library(latex2exp)
library(ggsci)
library(ggpubr)
output_dir <- ""
n <- 150
true_clusters <- c(rep(1, 50), rep(2, 50), rep(3, 50))
delta <- 0
q_seq <- c(5, 10, 25)
rho_seq <- c(0, 0.4, 0.8)

for (q in q_seq) {
  for (rho in rho_seq) {
    mu <- rbind(
      c(delta / 2, rep(0, q - 1)),
      c(rep(0, q - 1), delta / 2), c(rep(0, q - 1), delta / 2)
    )
    sig <- 1

    cov_mat_sim <- matrix(rho, nrow = q, ncol = q)
    diag(cov_mat_sim) <- 1
    cov_mat_sim <- cov_mat_sim * sig

    n_repeat <- 1500
    p_val_collection_hier <- c()
    p_val_collection_hier_centroid <- c()
    p_val_collection <- c()
    p_val_collection_hier_single <- c()
    p_naive_collection <- c()
    feat_collection <- c()

    for (i in c(1:n_repeat)) {
      cat("q=", q, "rho=", rho, "i=", i, "\n")
      set.seed(i)
      # exchangeable
      X <- mvrnorm(n, mu = rep(0, q), Sigma = cov_mat_sim) + mu[true_clusters, ]
      cluster_num <- sample(c(1:3), 3)
      feat <- sample(c(1:q), 1)
      cl_1_2_inference_demo <- kmeans_inference_1f(X,
        k = 3, cluster_num[1],
        cluster_num[2],
        feat = feat, iso = FALSE,
        sig = NULL, covMat = cov_mat_sim,
        iter.max = 15
      )
      p_val_collection[i] <- cl_1_2_inference_demo$pval
      p_naive_collection[i] <- cl_1_2_inference_demo$p_naive
      feat_collection[i] <- feat
      # single linkage is not gonna be great
      # centroid; average <- power and type i error
      hcl <- hclust(dist(X, method = "euclidean")^2, method = "average")
      result_hier_avg <- test_hier_clusters_exact_1f(
        X = X, link = "average", hcl = hcl, K = 3, k1 = cluster_num[1],
        k2 = cluster_num[2],
        indpt = F,
        covMat = cov_mat_sim,
        feat = feat
      )
      hcl <- hclust(dist(X, method = "euclidean")^2, method = "centroid")
      result_hier_centroid <- test_hier_clusters_exact_1f(
        X = X, link = "centroid", hcl = hcl, K = 3,
        k1 = cluster_num[1], k2 = cluster_num[2],
        indpt = F,
        covMat = cov_mat_sim,
        feat = feat
      )
      p_val_collection_hier[i] <- result_hier_avg$pval
      p_val_collection_hier_centroid[i] <- result_hier_centroid$pval

      hcl <- hclust(dist(X, method = "euclidean")^2, method = "single")
      result_hier_single <- test_hier_clusters_exact_1f(
        X = X, link = "single", hcl = hcl, K = 3,
        k1 = cluster_num[1], k2 = cluster_num[2],
        indpt = F,
        covMat = cov_mat_sim,
        feat = feat
      )
      p_val_collection_hier_single[i] <- result_hier_single$pval
    }

    if (delta > 0) {
      df_p0_plot <- data.frame(
        p_selective_hier = c(p_val_collection_hier),
        p_val_collection_hier_centroid = c(p_val_collection_hier_centroid),
        p_val_collection_hier_single = p_val_collection_hier_single,
        p_selective = c(p_val_collection),
        p_naive = c(p_naive_collection),
        feature_num = feat_collection,
        q = q,
        n = n
      ) %>%
        filter(feature_num != q & feature_num != 1) %>%
        dplyr::select(-c(feature_num, q, n)) %>%
        pivot_longer(everything(), names_to = "p_type", values_to = "p_values") %>%
        mutate(p_type = as.factor(p_type)) %>%
        group_by(p_type) %>%
        mutate(theoretical = ecdf(p_values)(p_values)) %>%
        mutate(q = q, n = n, rho = rho)
    } else {
      df_p0_plot <- data.frame(
        p_selective_hier = c(p_val_collection_hier),
        p_val_collection_hier_centroid = c(p_val_collection_hier_centroid),
        p_selective = c(p_val_collection),
        p_val_collection_hier_single = p_val_collection_hier_single,
        p_naive = c(p_naive_collection),
        feature_num = feat_collection,
        q = q,
        n = n
      ) %>%
        dplyr::select(-c(feature_num, q, n)) %>%
        pivot_longer(everything(), names_to = "p_type", values_to = "p_values") %>%
        mutate(p_type = as.factor(p_type)) %>%
        group_by(p_type) %>%
        mutate(theoretical = ecdf(p_values)(p_values)) %>%
        mutate(q = q, n = n, rho = rho)
    }


    save(df_p0_plot, file = paste0(output_dir, "single_feature_type_1_q_", q, "_n_", n, "_rho_", rho, ".RData"))
  }
}
