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
  pattern = glob2rx("power_single_feature_type_1_q_*_rho_*.*"), full.names = T
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

rho_label_names <- c(
  `0` = ("cor=0"),
  `0.4` = ("cor=0.4"),
  `0.8` = ("cor=0.8")
)




detect_prob_df <- df_p_value_list %>%
  filter((delta %in% c(3, 4, 5, 6, 7, 8)) & (q %in% c(5, 10, 25))) %>%
  group_by(q, rho, delta) %>%
  summarise(
    mean_detect_p_k = mean(aRI_seq == 1, na.rm = T),
    sd_detect_p_k = sqrt(mean_detect_p_k * (1 - mean_detect_p_k) / n()),
    mean_detect_p_avg = mean(aRI_seq_avg == 1, na.rm = T),
    sd_detect_p_avg = sqrt(mean_detect_p_avg * (1 - mean_detect_p_avg) / n()),
    mean_detect_p_centroid = mean(aRI_seq_centroid == 1, na.rm = T),
    sd_detect_p_centroid = sqrt(mean_detect_p_centroid * (1 - mean_detect_p_centroid) / n()),
    mean_detect_p_single = mean(aRI_seq_single == 1, na.rm = T),
    sd_detect_p_single = sqrt(mean_detect_p_single * (1 - mean_detect_p_single) / n())
  )


# Reshape the dataframe to have separate columns for means and standard deviations
detect_prob_df_flat <- detect_prob_df %>%
  pivot_longer(
    cols = starts_with(c("mean_", "sd_")),
    names_to = c(".value", "method"),
    names_pattern = "(mean|sd)_(.*)"
  )

# Rename the 'method' column values for better clarity
detect_prob_df_flat$method <- sub("detect_p_", "", detect_prob_df_flat$method)
appender <- function(string) {
  TeX(paste(r'($\lambda$)', string))
}

# rho_label_names <- c(
#  `0` = 'cor=0',
#  `0.4` = 'cor=0.4',
#  `0.8` = 'cor=0.8')

# try parsing as well
detect_prob_df_flat$rho <- factor(detect_prob_df_flat$rho, labels = c(
  "0" = parse(text = TeX("$\\rho$=0")),
  "0.4" = parse(text = TeX("$\\rho$=0.4")),
  "0.8" = parse(text = TeX("$\\rho$=0.8"))
))

q_label_names <- c(
  `5` = "q=5",
  `10` = "q=10",
  `25` = "q=25"
)

# levels(detect_prob_df_flat$rho_annotation) <- c("0" = latex2exp::TeX("$\\rho$=0",output = 'expression'),
#                                     "0.4" = latex2exp::TeX("$\\rho$=0.4",output = 'character'),
#                                     "0.8" = latex2exp::TeX("$\\rho$=0.8",output = 'character'))

p_detect_prob <- detect_prob_df_flat %>%
  ggplot(aes(
    x = delta, y = mean, color = as.factor(method),
    # linetype=as.factor(rho),
    group = interaction(method, rho)
  )) +
  geom_point() +
  geom_pointrange(aes(
    ymin = mean + sd,
    ymax = mean - sd
  )) +
  geom_line() +
  facet_grid(q ~ rho,
    scales = "free",
    labeller = labeller(.rows = q_label_names, .cols = label_parsed, )
  ) +
  ylab(TeX("Detection probability")) +
  xlab(TeX("Distance between features in true clusters: $\\delta$")) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  labs(color = "") +
  scale_color_manual(
    values = c("#00A1D5FF", "#B24745FF", "#DF8F44FF", "#79AF97FF"),
    labels = unname(TeX(c(
      "Average linkage",
      "Centroid linkage",
      "K-means",
      "Single linkage"
    )))
  )


png(paste0("~/Desktop/fall_2022/p_detect_prob.png"),
  width = 12, height = 11, res = 400, units = "in"
)
p_detect_prob
dev.off()


p_detect_prob_aistat <- detect_prob_df_flat %>%
  filter(q == 10) %>%
  ggplot(aes(
    x = delta, y = mean, color = as.factor(method),
    # linetype=as.factor(rho),
    group = interaction(method, rho)
  )) +
  geom_point() +
  geom_pointrange(aes(
    ymin = mean + sd,
    ymax = mean - sd
  )) +
  geom_line() +
  facet_grid(~rho,
    scales = "free",
    labeller = labeller(.cols = label_parsed, )
  ) +
  ylab(TeX("Detection probability")) +
  xlab(TeX("Distance between features in true clusters: $\\delta$")) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  labs(color = "") +
  scale_color_manual(
    values = c("#00A1D5FF", "#B24745FF", "#DF8F44FF", "#79AF97FF"),
    labels = unname(TeX(c(
      "Average linkage",
      "Centroid linkage",
      "K-means",
      "Single linkage"
    )))
  )


png(paste0("~/Desktop/fall_2022/p_detect_prob.png"),
  width = 12, height = 11, res = 400, units = "in"
)
p_detect_prob
dev.off()



alpha_thresh <- 0.05


colnames(df_p_value_list) <- c(
  "p_avg", "p_centroid", "p_kmeans", "feat_kmeans",
  "feat_centroid", "feat_avg", "aRI_centroid", "aRI_avg",
  "aRI_kmeans", "q", "n", "rho", "delta", "p_single",
  "feat_single", "aRI_single"
)

df_p_value_list_long <- df_p_value_list %>%
  filter((delta %in% c(3, 4, 5, 6, 7, 8)) & (q %in% c(5, 10, 25))) %>%
  pivot_longer(
    cols = starts_with(c("p_", "feat_", "aRI_")),
    names_to = c(".value", "method"),
    names_pattern = "(p|feat|aRI)_(.*)"
  )

# we need to fix the plot -- double feature being counted right now
cond_power_df <- df_p_value_list_long %>%
  # mutate(effect_size = abs(feat)) %>%
  # filter(effect_size==delta/2|effect_size==delta) %>%
  group_by(q, rho, delta, method) %>%
  summarise(
    mean_power = sum((p <= alpha_thresh) & (aRI == 1)) / sum(aRI == 1),
    sd_power = sqrt(mean_power * (1 - mean_power) / sum(aRI == 1))
  ) %>%
  mutate(
    across(c(mean_power, sd_power), ~ replace_na(.x, 0))
  )

df_p_value_list_long$rho <- factor(df_p_value_list_long$rho, labels = c(
  "0" = parse(text = TeX("$\\rho$=0")),
  "0.4" = parse(text = TeX("$\\rho$=0.4")),
  "0.8" = parse(text = TeX("$\\rho$=0.8"))
))

smooth_power <- df_p_value_list_long %>%
  mutate(effect_size = abs(feat)) %>%
  ggplot(aes(
    x = effect_size,
    y = as.numeric(p <= alpha_thresh),
    colour = method
  )) +
  geom_smooth() +
  facet_grid(q ~ rho,
    scales = "free",
    labeller = labeller(.rows = q_label_names, .cols = label_parsed, )
  ) +
  ylab(TeX("Power at $\\alpha = 0.05$")) +
  xlab(TeX("$ [\\mu^T\\nu]_j $")) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  labs(color = "Test") +
  scale_color_manual(
    values = c("#00A1D5FF", "#B24745FF", "#DF8F44FF", "#79AF97FF"),
    labels = unname(TeX(c(
      "$p_{j,average}$", "$p_{j,centroid}$",
      "$p_{j,kmeans}$", "$p_{j,single}"
    )))
  )




p_smooth_power_aistat <- df_p_value_list_long %>%
  mutate(effect_size = abs(feat)) %>%
  filter(q == 10) %>%
  ggplot(aes(
    x = effect_size,
    y = as.numeric(p <= alpha_thresh),
    colour = method
  )) +
  geom_smooth() +
  facet_grid(. ~ rho,
    scales = "free",
    labeller = labeller(.cols = label_parsed)
  ) +
  ylab(TeX("Power at $\\alpha = 0.05$")) +
  xlab(TeX("$ [\\mu^T\\nu]_j $")) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  labs(color = "Test") +
  scale_color_manual(
    values = c("#00A1D5FF", "#B24745FF", "#DF8F44FF", "#79AF97FF"),
    labels = unname(TeX(c(
      "$p_{j,average}$", "$p_{j,centroid}$",
      "$p_{j,kmeans}$", "$p_{j,single}"
    )))
  )


png(paste0("~/Desktop/fall_2022/p_smooth_power_aistat.png"),
  width = 8, height = 3.5, res = 400, units = "in"
)
p_smooth_power_aistat
dev.off()


png(paste0("~/Desktop/fall_2022/p_detect_prob_aistat.png"),
  width = 8, height = 3.5, res = 400, units = "in"
)
p_detect_prob_aistat
dev.off()



png(paste0("~/Desktop/fall_2022/Smooth_Power.png"),
  width = 12, height = 11, res = 400, units = "in"
)
smooth_power
dev.off()

cond_power_df$rho <- factor(cond_power_df$rho, labels = c(
  "0" = parse(text = TeX("$\\rho$=0")),
  "0.4" = parse(text = TeX("$\\rho$=0.4")),
  "0.8" = parse(text = TeX("$\\rho$=0.8"))
))


p_cond_power <- cond_power_df %>%
  ggplot(aes(x = delta, y = mean_power, color = as.factor(method), )) +
  geom_point() +
  geom_pointrange(aes(
    ymin = mean_power + sd_power,
    ymax = mean_power - sd_power
  )) +
  geom_line() +
  facet_grid(q ~ rho,
    scales = "free",
    labeller = labeller(.rows = q_label_names, .cols = label_parsed, )
  ) +
  ylab(TeX("Conditional power at $\\alpha = 0.05$")) +
  xlab(TeX("Distance between features in true clusters: $\\delta$")) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  labs(color = "Test") +
  scale_color_manual(
    values = c("#00A1D5FF", "#B24745FF", "#DF8F44FF", "#79AF97FF"),
    labels = unname(TeX(c(
      "$p_{j,average}$", "$p_{j,centroid}$",
      "$p_{j,kmeans}$", "$p_{j,single}"
    )))
  )



png(paste0("~/Desktop/fall_2022/Conditional_Power.png"),
  width = 12, height = 11, res = 400, units = "in"
)
p_cond_power
dev.off()

p_cond_power_aistat <- cond_power_df %>%
  filter(q == 10) %>%
  ggplot(aes(x = delta, y = mean_power, color = as.factor(method), )) +
  geom_point() +
  geom_pointrange(aes(
    ymin = mean_power + sd_power,
    ymax = mean_power - sd_power
  )) +
  geom_line() +
  facet_grid(. ~ rho,
    scales = "free",
    labeller = label_parsed
  ) +
  ylab(TeX("Conditional power at $\\alpha = 0.05$")) +
  xlab(TeX("Distance between features in true clusters: $\\delta$")) +
  theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "lines"),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15, hjust = 0),
    axis.title = element_text(size = 15)
  ) +
  labs(color = "Test") +
  scale_color_manual(
    values = c("#00A1D5FF", "#B24745FF", "#DF8F44FF", "#79AF97FF"),
    labels = unname(TeX(c(
      "$p_{j,average}$", "$p_{j,centroid}$",
      "$p_{j,kmeans}$", "$p_{j,single}"
    )))
  )



png(paste0("./Figure_4.png"),
  width = 8, height = 3.7, res = 400, units = "in"
)
p_cond_power_aistat
dev.off()
