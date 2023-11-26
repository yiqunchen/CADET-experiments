packages <- c(
  "devtools", "tidyverse", "latex2exp", "RColorBrewer",
  "shades", "patchwork", "ggsci", "ggpubr", "corrplot",
  "umap", "DropletUtils", "scater", "ggfortify", "mclust",
  "rdist", "MASS"
)
install.packages(setdiff(packages, rownames(installed.packages())))
if (!("CADET" %in% rownames(installed.packages()))) {
  devtools::install_github("yiqunchen/CADET")
}
