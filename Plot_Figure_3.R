set.seed(2022)
library(CADET)
library(MASS)
library(corrplot)
library(tidyverse)
library(latex2exp)
library(ggsci)
library(ggpubr)

input_dir <- ""
sim_files <- list.files(
  path = input_dir,
  pattern = glob2rx("single_feature_type_1_q_*_rho_*.*"), full.names = T
)


concat_p_value_list <- list()

for (i in seq_along(sim_files)) {
  current_file <- sim_files[i]
  load(file = current_file)
  concat_p_value_list[[i]] <- df_p0_plot
  rm(df_p0_plot)
}

df_p_value_list <- dplyr::bind_rows(concat_p_value_list)

alpha_thresh <- 0.05
# load power data


q_label_names <- c(
  `5` = "q=5",
  `10` = "q=10",
  `25` = "q=25"
)


# try parsing as well
df_p_value_list$rho <- factor(df_p_value_list$rho, labels = c(
  "0" = parse(text = TeX("$\\rho$=0")),
  "0.4" = parse(text = TeX("$\\rho$=0.4")),
  "0.8" = parse(text = TeX("$\\rho$=0.8"))
))


p_aggregate_type_1 <- df_p_value_list %>%
  filter(q %in% c(5, 10, 25)) %>%
  # mutate(p_type = fct_relevel(as.factor(p_type), "p_naive", after = Inf)) %>%
  group_by(q, p_type, rho) %>%
  mutate(theoretical = ecdf(p_values)(p_values)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(y = p_values, x = (theoretical), colour = p_type), size = 0.6) +
  facet_grid(q ~ rho, scales = "free", labeller = labeller(q = q_label_names, rho = label_parsed)) +
  ylab("P-value quantiles") +
  xlab("Uniform(0,1) quantiles") +
  # ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(color = "Test") +
  scale_color_jama(labels = unname(TeX(c(
    "$p_{j,naive}$", "$p_{j,kmeans}$",
    "$p_{j,average}$", "$p_{j,centroid}$"
  )))) +
  theme_bw(base_size = 17) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    legend.position = "bottom",
    legend.title = element_text(size = 17),
    axis.text = element_text(size = 13),
    legend.text = element_text(size = 17, hjust = 0),
    axis.title = element_text(size = 17)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5)))


png(paste0("~/Desktop/fall_2022/Type_1_error_control.png"),
  width = 12, height = 12.5, res = 400, units = "in"
)
p_aggregate_type_1
dev.off()

p_aistat_paper_main <- df_p_value_list %>%
  filter(q %in% c(10)) %>%
  # mutate(p_type = fct_relevel(as.factor(p_type), "p_naive", after = Inf)) %>%
  group_by(q, p_type, rho) %>%
  mutate(theoretical = ecdf(p_values)(p_values)) %>%
  ungroup() %>%
  ggplot() +
  geom_point(aes(y = p_values, x = (theoretical), colour = p_type), size = 0.6) +
  facet_grid(. ~ rho, scales = "free", labeller = label_parsed) +
  ylab("P-value quantiles") +
  xlab("Uniform(0,1) quantiles") +
  # ggtitle(TeX('Selective test p-values'))+
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(color = "Test") +
  scale_color_jama(labels = unname(TeX(c(
    "$p_{j,naive}$", "$p_{j,kmeans}$",
    "$p_{j,average}$", "$p_{j,centroid}$",
    "$p_{j,single}$"
  )))) +
  theme_bw(base_size = 17) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 17),
    legend.position = "bottom",
    legend.title = element_text(size = 17),
    axis.text = element_text(size = 13),
    legend.text = element_text(size = 17, hjust = 0),
    axis.title = element_text(size = 17)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5)))


png(paste0("./Figure_3.png"),
  width = 8, height = 3.5, res = 400, units = "in"
)
p_aistat_paper_main
dev.off()
