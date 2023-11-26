set.seed(2022)
library(CADET)
library(MASS)
library(corrplot)
library(tidyverse)
library(latex2exp)
library(ggpubr)
library(mclust)

output_dir <- ""
n <- 150
n_repeat <- 2000
true_clusters <- c(rep(1, 50), rep(2, 50), rep(3, 50))
q_seq <- c(5, 10, 25)
rho_seq <- c(0, 0.4, 0.8)
delta_seq <- c(3, 4, 5, 6, 7, 8)
for (delta in delta_seq) {
  for (q in q_seq) {
    for (rho in rho_seq) {
      mu <- rbind(
        c(rep(0, floor(q / 2)), rep(delta / 2, q - floor(q / 2))),
        rep(0, q),
        c(rep(0, floor(q / 2)), rep(-delta / 2, q - floor(q / 2)))
      )
      sig <- 1

      cov_mat_sim <- matrix(rho, nrow = q, ncol = q)
      diag(cov_mat_sim) <- 1
      cov_mat_sim <- cov_mat_sim * sig
      p_val_collection_hier <- c()
      p_val_collection_hier_centroid <- c()
      p_val_collection_hier_single <- c()
      p_val_collection <- c()
      p_naive_collection <- c()
      feat_collection_kmeans <- c()
      feat_collection_hier_avg <- c()
      feat_collection_hier_centroid <- c()
      feat_collection_hier_single <- c()
      aRI_seq <- c()
      aRI_seq_avg <- c()
      aRI_seq_centroid <- c()
      aRI_seq_single <- c()


      for (i in c(1:n_repeat)) {
        if (i %% 50 == 0) {
          cat("delta=", delta, "q=", q, "rho=", rho, "i=", i, "\n")
        }
        set.seed(i)
        # exchangeable
        mu_mat <- mu[true_clusters, ]
        X <- mvrnorm(n, mu = rep(0, q), Sigma = cov_mat_sim) + mu_mat
        cluster_num <- sample(c(1:3), 3)
        feat <- sample(c((floor(q / 2) + 1):q), 1)
        cl_1_2_inference_demo <- kmeans_inference_1f(X,
          k = 3, cluster_num[1],
          cluster_num[2],
          feat = feat, iso = FALSE,
          sig = NULL, covMat = cov_mat_sim,
          iter.max = 15
        )
        p_val_collection[i] <- cl_1_2_inference_demo$pval
        final_cluster_kmeans <- cl_1_2_inference_demo$final_cluster
        aRI_seq[i] <- mclust::adjustedRandIndex(final_cluster_kmeans, true_clusters)
        feat_collection_kmeans[i] <- mean(mu_mat[final_cluster_kmeans == cluster_num[1], feat]) - mean(mu_mat[final_cluster_kmeans == cluster_num[2], feat])
        # single linkage is not gonna be great
        # centroid; average <- power and type i error
        hcl <- hclust(dist(X, method = "euclidean")^2, method = "average")
        cut_avg <- cutree(hcl, k = 3)
        aRI_seq_avg[i] <- mclust::adjustedRandIndex(cut_avg, true_clusters)
        result_hier_avg <- test_hier_clusters_exact_1f(
          X = X, link = "average", hcl = hcl, K = 3, k1 = cluster_num[1],
          k2 = cluster_num[2],
          indpt = F,
          covMat = cov_mat_sim,
          feat = feat
        )
        feat_collection_hier_avg[i] <- mean(mu_mat[cut_avg == cluster_num[1], feat]) - mean(mu_mat[cut_avg == cluster_num[2], feat])

        hcl <- hclust(dist(X, method = "euclidean")^2, method = "centroid")


        result_hier_centroid <- test_hier_clusters_exact_1f(
          X = X, link = "centroid", hcl = hcl, K = 3,
          k1 = cluster_num[1], k2 = cluster_num[2],
          indpt = F,
          covMat = cov_mat_sim,
          feat = feat
        )
        cut_centroid <- cutree(hcl, k = 3)
        aRI_seq_centroid[i] <- mclust::adjustedRandIndex(cut_centroid, true_clusters)
        feat_collection_hier_centroid[i] <- mean(mu_mat[cut_centroid == cluster_num[1], feat]) - mean(mu_mat[cut_centroid == cluster_num[2], feat])
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
        cut_centroid <- cutree(hcl, k = 3)
        aRI_seq_single[i] <- mclust::adjustedRandIndex(cut_centroid, true_clusters)
        feat_collection_hier_single[i] <- mean(mu_mat[cut_centroid == cluster_num[1], feat]) - mean(mu_mat[cut_centroid == cluster_num[2], feat])
        p_val_collection_hier_single[i] <- result_hier_single$pval
      }



      df_p0_plot <- data.frame(
        p_selective_hier = c(p_val_collection_hier),
        p_val_collection_hier_centroid = c(p_val_collection_hier_centroid),
        p_selective = c(p_val_collection),
        p_val_collection_hier_single = p_val_collection_hier_single,
        feat_collection_kmeans = feat_collection_kmeans,
        feat_collection_hier_centroid = feat_collection_hier_centroid,
        feat_collection_hier_avg = feat_collection_hier_avg,
        feat_collection_hier_single = feat_collection_hier_single,
        aRI_seq_centroid = aRI_seq_centroid,
        aRI_seq_avg = aRI_seq_avg,
        aRI_seq_single = aRI_seq_single,
        aRI_seq = aRI_seq
      ) %>%
        mutate(q = q, n = n, rho = rho, delta = delta)

      save(df_p0_plot, file = paste0(output_dir, "power_single_feature_type_1_q_", q, "_n_", n, "_rho_", rho, "_delta_", delta, ".RData"))
    }
  }
}
