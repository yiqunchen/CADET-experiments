library(Seurat)
library(DropletUtils)
library(scater)
library(ggfortify)
library(patchwork)
library(CADET)
library(umap)
library(SeuratDisk)
library(zellkonverter)
library(scuttle)
### qq plot
library(tidyverse)
library(latex2exp)
library(ggsci)
library(ggpubr)
library(CADET)

set.seed(2022)

TS_bone_marrow <- zellkonverter::readH5AD("~/Downloads/TS_Bone_Marrow.h5ad",
  verbose = T,
  layers = F, varm = F, obsm = F, varp = F, obsp = F, uns = F
)
# all manually annotated
table(TS_bone_marrow$manually_annotated)
names(assays(TS_bone_marrow)) <- "counts"

process_data <- function(sce) {
  # Define mitochondrial, ribosomal protein genes
  rowData(sce)$Mito <- grepl("^MT", rowData(sce)$Symbol)
  rowData(sce)$ProteinRibo <- grepl("^RP", rowData(sce)$Symbol)

  # Calculate quality statistics
  sce <- addPerCellQC(sce, subsets = list(mt = rowData(sce)$Mito, rbp = rowData(sce)$ProteinRibo))
  sce <- addPerFeatureQC(sce)

  # Delete bad cells (low library size, low gene count, high mitochondrial read proportion)
  libsize.drop <- isOutlier(sce$sum, nmads = 3, type = "lower", log = TRUE)
  feature.drop <- isOutlier(sce$detected, nmads = 3, type = "both", log = TRUE)
  mito.drop <- isOutlier(sce$subsets_mt_percent, nmads = 3, type = "higher", log = TRUE)
  sce <- sce[, !(libsize.drop | feature.drop | mito.drop)]

  # Normalize library sizes, then log2 transformation with pseudo-count of 1
  sce <- logNormCounts(sce)

  # Subset to top 500 average count genes
  X <- t(logcounts(sce))
  X <- as.matrix(X)
  X <- X[, order(rowData(sce)$mean, decreasing = T)[1:500]]
  return(list(data = X, labels = colData(sce)$CellType))
}

X_bone_marrow <- t(logcounts(logNormCounts(TS_bone_marrow)))
organ_tissue <- TS_bone_marrow$organ_tissue

X1 <- X_bone_marrow[which(
  (TS_bone_marrow$donor == "TSP14") &
    (TS_bone_marrow$cell_ontology_class == "cd4-positive, alpha-beta t cell")
), ]
# throw an error for zero var genes
var_X1 <- apply(X1, 2, var)
mean_X1 <- apply(X1, 2, mean)
X1 <- X1[, order(var_X1, decreasing = T)[1:500]]
# prop of zeros
zero_prop_X1 <- apply(X1, 2, mean)
X1 <- as.matrix(X1)
summary(apply(X1, 2, var))
hist(c(X1))
# X1 <- X.endo.ss
cov_sample <- coop::covar(X1)
set.seed(1234)
seed_list <- sample(2021, size = 30, replace = FALSE)
k_means_list <- lapply(seed_list, function(x) {
  kmeans_estimation(X1,
    k = 4, 1, iter.max = 30, seed = x
  )
})
within_ss_list <- lapply(k_means_list, function(x) x$objective[[x$iter]])
best_seed <- seed_list[which.min(unlist(within_ss_list))]
current_kmeans_k_3_neg <- k_means_list[[which.min(unlist(within_ss_list))]]
final_cluster_k_3_neg <- current_kmeans_k_3_neg$cluster[[current_kmeans_k_3_neg$iter]]
table(final_cluster_k_3_neg)

set.seed(2021)
X1_umap <- umap(X1)
umap_data_x1 <- X1_umap$layout
colnames(umap_data_x1) <- c("UMAP1", "UMAP2")


# accumulate all the p-values
selective_pvals_c1_c2 <- rep(NA, dim(X1)[2])
selective_pvals_c1_c3 <- rep(NA, dim(X1)[2])
selective_pvals_c2_c3 <- rep(NA, dim(X1)[2])
selective_pvals_c2_c4 <- rep(NA, dim(X1)[2])
selective_pvals_c1_c4 <- rep(NA, dim(X1)[2])
selective_pvals_c3_c4 <- rep(NA, dim(X1)[2])

naive_pvals_c1_c2 <- rep(NA, dim(X1)[2])
naive_pvals_c1_c3 <- rep(NA, dim(X1)[2])
naive_pvals_c2_c3 <- rep(NA, dim(X1)[2])
naive_pvals_c2_c4 <- rep(NA, dim(X1)[2])
naive_pvals_c1_c4 <- rep(NA, dim(X1)[2])
naive_pvals_c3_c4 <- rep(NA, dim(X1)[2])

for (i in c(1:dim(X1)[2])) {
  cat("i", i, "\n")

  cl_1_4_feat_i <- kmeans_inference_1f(X1,
    k = 4, 1, 4,
    feat = i, iso = FALSE, sig = NULL,
    covMat = cov_sample, seed = best_seed,
    iter.max = 30
  )
  selective_pvals_c1_c4[i] <- cl_1_4_feat_i$pval
  naive_pvals_c1_c4[i] <- cl_1_4_feat_i$p_naive

  cl_3_4_feat_i <- kmeans_inference_1f(X1,
    k = 4, 3, 4,
    feat = i, iso = FALSE, sig = NULL,
    covMat = cov_sample, seed = best_seed,
    iter.max = 30
  )
  selective_pvals_c3_c4[i] <- cl_3_4_feat_i$pval
  naive_pvals_c3_c4[i] <- cl_3_4_feat_i$p_naive

  cl_1_2_feat_i <- kmeans_inference_1f(X1,
    k = 4, 1, 2,
    feat = i, iso = FALSE, sig = NULL,
    covMat = cov_sample, seed = best_seed,
    iter.max = 30
  )
  selective_pvals_c1_c2[i] <- cl_1_2_feat_i$pval
  naive_pvals_c1_c2[i] <- cl_1_2_feat_i$p_naive

  cl_2_3_feat_i <- kmeans_inference_1f(X1,
    k = 4, 2, 3,
    feat = i, iso = FALSE, sig = NULL,
    covMat = cov_sample, seed = best_seed,
    iter.max = 30
  )
  selective_pvals_c2_c3[i] <- cl_2_3_feat_i$pval
  naive_pvals_c2_c3[i] <- cl_2_3_feat_i$p_naive


  cl_1_3_feat_i <- kmeans_inference_1f(X1,
    k = 4, 1, 3,
    feat = i, iso = FALSE, sig = NULL,
    covMat = cov_sample, seed = best_seed,
    iter.max = 30
  )
  selective_pvals_c1_c3[i] <- cl_1_3_feat_i$pval
  naive_pvals_c1_c3[i] <- cl_1_3_feat_i$p_naive

  cl_2_4_feat_i <- kmeans_inference_1f(X1,
    k = 4, 2, 4,
    feat = i, iso = FALSE, sig = NULL,
    covMat = cov_sample, seed = best_seed,
    iter.max = 30
  )
  selective_pvals_c2_c4[i] <- cl_2_4_feat_i$pval
  naive_pvals_c2_c4[i] <- cl_2_4_feat_i$p_naive
}


plot_output_dir <- ""



df_p0_summary <- data.frame(
  p_select = c(
    selective_pvals_c1_c2, selective_pvals_c2_c3,
    selective_pvals_c1_c3,
    selective_pvals_c2_c4,
    selective_pvals_c3_c4, selective_pvals_c1_c4
  ),
  p_naive = c(
    naive_pvals_c1_c2, naive_pvals_c2_c3,
    naive_pvals_c1_c3, naive_pvals_c2_c4,
    naive_pvals_c1_c4, naive_pvals_c3_c4
  ),
  cluster_1 = rep(c(1, 2, 1, 2, 3, 1), each = dim(X1)[2]),
  cluster_2 = rep(c(2, 3, 3, 4, 4, 4), each = dim(X1)[2])
) %>%
  mutate(cluster_pair = paste0(cluster_1, "-", cluster_2))

df_p0_plot <- df_p0_summary %>%
  pivot_longer(-c(cluster_1, cluster_2, cluster_pair), names_to = "p_type", values_to = "p_values") %>%
  mutate(p_type = as.factor(p_type)) %>%
  group_by(p_type) %>%
  mutate(theoretical = ecdf(p_values)(p_values))

type_1_plot <- ggplot(data = df_p0_plot) +
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
    legend.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 23, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  labs(color = "Test") +
  scale_color_manual(
    values = c("#374E55FF", "#DF8F44FF"),
    labels = unname(TeX(c("$\\hat{p}_{j,naive}$", "$\\hat{p}_{j,k-means}$")))
  )


png("Figure_6_a.png",
  width = 7.5, height = 7.5, res = 300, units = "in"
)
type_1_plot
dev.off()

# try to understand which genes are differentiually expressed
# df_p0_summary$gene_name <- rep(colnames(X1),6)
group_BH <- df_p0_summary %>%
  group_by(cluster_pair) %>%
  group_modify(~ broom::tidy(p.adjust(.x$p_select, method = "BH"))) %>%
  pull(x)

group_BH_naive <- df_p0_summary %>%
  group_by(cluster_pair) %>%
  group_modify(~ broom::tidy(p.adjust(.x$p_naive, method = "BH"))) %>%
  pull(x)
# try to understand which genes are differentiually expressed
bh_adjusted_select_pvalues <- p.adjust(df_p0_summary$p_select, method = "BH")
bh_adjusted_naive_pvalues <- p.adjust(df_p0_summary$p_naive, method = "BH")
## plot BH against FDR control
df_p0_summary$bh_adjusted_select_pvalues_by_pair <- group_BH
df_p0_summary$bh_adjusted_naive_pvalues_by_pair <- group_BH_naive
fdr_cutoff <- seq(0, 1, by = 0.05)
cluster_pairs <- unique(df_p0_summary$cluster_pair)
# agg data by cluster pair
plot_fdr_null_num_rej_by_pair <- data.frame()
for (cluster in cluster_pairs) {
  fdr_selective_rej <- sapply(fdr_cutoff, function(x) {
    sum(df_p0_summary %>% filter(cluster_pair == cluster) %>% pull(bh_adjusted_select_pvalues_by_pair) <= x)
  })
  fdr_naive_rej <- sapply(fdr_cutoff, function(x) {
    sum(df_p0_summary %>% filter(cluster_pair == cluster) %>% pull(bh_adjusted_naive_pvalues_by_pair) <= x)
  })
  plot_fdr_null_num_rej <- data.frame(
    fdr_cutoff, fdr_selective_rej, fdr_naive_rej,
    cluster
  )
  plot_fdr_null_num_rej_by_pair <- rbind(plot_fdr_null_num_rej_by_pair, plot_fdr_null_num_rej)
}

plot_fdr_null_num_rej_by_pair_df <- plot_fdr_null_num_rej_by_pair %>%
  pivot_longer(c(fdr_selective_rej, fdr_naive_rej),
    names_to = "method",
    values_to = "num_of_rejections"
  )

plot_fdr_null_num_rej_df <- plot_fdr_null_num_rej %>%
  pivot_longer(-c(fdr_cutoff, cluster),
    names_to = "method",
    values_to = "num_of_rejections"
  )


png("Figure_6_b.png",
  width = 12, height = 8, res = 300, units = "in"
)
ggplot(plot_fdr_null_num_rej_by_pair_df %>% filter(fdr_cutoff <= 0.25), aes(
  x = fdr_cutoff,
  y = num_of_rejections, colour = method
), ) +
  geom_line(size = 1.5) +
  ylab("Number of rejected null hypotheses") +
  xlab("Nominal FDR levels") +
  # ggtitle("Bone marrow CD 4 alpha/beta, Donor 14")+
  labs(color = "Test") +
  scale_color_jama(labels = unname(TeX(c(
    "BH-adjusted $\\hat{p}_{j,naive}$",
    "BH-ajusted $\\hat{p}_{j,k-means}$"
  )))) +
  scale_linetype_manual(values = c("solid", "solid")) +
  # scale_size_discrete(position="none")+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    axis.text = element_text(size = 9),
    legend.text = element_text(size = 20, hjust = 0),
    axis.title = element_text(size = 12)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  facet_wrap(. ~ cluster, nrow = 2)
dev.off()
