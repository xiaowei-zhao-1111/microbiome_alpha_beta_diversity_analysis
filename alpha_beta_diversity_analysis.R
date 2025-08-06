### Load required libraries
library(lmerTest)      # For mixed-effects models
library(vegan)         # For diversity indices and PERMANOVA
library(ade4)          # For PCoA
library(ggplot2)       # For visualization
library(reshape2)      # For data reshaping
library(dplyr)         # For data manipulation
library(ggpubr)        # For publication-ready plots
library(tidyr)         # For tidying data
library(RColorBrewer)  # For color palettes

###############################################################################
# Set project folder and load data
folder_path <- "~/alpha_beta_diversity_project"
setwd(folder_path)

# Read species relative abundance table from MetaPhlAn output
species_data <- read.table(file = "~/alpha_beta_diversity_project/species_abundance.tsv", sep = "\t", header = TRUE)
species <- species_data[grepl("s__", species_data$clade_name) & !grepl("t__", species_data$clade_name), ]

# Clean species names
rownames(species) <- gsub(".*s__", "", species$clade_name)
species <- species[, -1]

# Standardize sample names
colnames(species) <- gsub("metaphlan_|_S.*", "", names(species))

# Load metadata
metadata <- read.csv(file = "~/alpha_beta_diversity_project/metadata.csv", sep = ",", header = TRUE, row.names = 1)
species_clean <- species[, rownames(metadata)]

# Convert relative abundances to proportions
species_prop <- data.frame(apply(species_clean, 2, function(x) x / sum(x)))

# Determine cutoff for low-abundance species using rank plot
relab_number <- unlist(species_prop, use.name = FALSE)
relab_number_ordered <- relab_number[order(relab_number, decreasing = TRUE)]

pdf("rank_plot_species_cutoff.pdf", height = 6, width = 8)
plot(x = 1:length(relab_number_ordered), y = log10(relab_number_ordered), pch = 20,
     xlab = "Rank", ylab = "Relative abundance", main = "Rank plot", xlim = c(0, 14000))
abline(h = -5, col = "red")
abline(h = -4, col = "blue")
abline(h = -4.5, col = "orange")
dev.off()

cut_off <- 10^-4.5
species_cutoff <- species_prop
species_cutoff[species_cutoff <= cut_off] <- 0

# Merge metadata into abundance table
species_t <- as.data.frame(t(species_cutoff))
species_all <- merge(species_t, metadata, by = "row.names", all = FALSE)
rownames(species_all) <- species_all$Row.names
species_all <- species_all[, -1]

###############################################################################
## Alpha- and Beta-diversity analysis

# Create output directory
folder_name <- "alpha_diversity_results"
folder <- file.path(folder_path, folder_name)
if (!file.exists(folder)) dir.create(folder)

# Subset to case-control comparisons
groupA_vs_control_data <- species_all[species_all$condition %in% c("groupA", "control"), ]

# Compute alpha diversity metrics
shannon <- diversity(groupA_vs_control_data[, 1:(ncol(groupA_vs_control_data)-32)], index = "shannon")
simpson <- diversity(groupA_vs_control_data[, 1:(ncol(groupA_vs_control_data)-32)], index = "simpson")
invsimpson <- diversity(groupA_vs_control_data[, 1:(ncol(groupA_vs_control_data)-32)], index = "invsimpson")
richness <- apply(groupA_vs_control_data[, 1:(ncol(groupA_vs_control_data)-32)], 1, function(x) sum(x > 0))

div_df <- data.frame(
  condition = groupA_vs_control_data$condition,
  household_id = groupA_vs_control_data$household_id,
  age = groupA_vs_control_data$age,
  BMI = groupA_vs_control_data$BMI,
  sex = groupA_vs_control_data$sex,
  richness = richness,
  shannon = shannon,
  simpson = simpson,
  invsimpson = invsimpson
)

write.csv(div_df, file.path(folder, "alpha_diversity_groupA_vs_control.csv"), row.names = TRUE)

# Mixed-effects models (adjusting for household)
diversity_list <- c("shannon", "simpson", "invsimpson", "richness")
p_values <- sapply(diversity_list, function(metric) {
  summary(lmer(as.formula(paste(metric, "~ condition + (1|household_id)")), data = div_df, REML = FALSE))["coefficients"][[1]][2,5]
})
mixed_effect_p <- data.frame(diversity = diversity_list, p_value = p_values)
write.csv(mixed_effect_p, file.path(folder, "mixed_effect_p_groupA_vs_control.csv"))

# Wilcoxon tests
wilcoxon_p <- sapply(diversity_list, function(metric) {
  with(div_df, wilcox.test(get(metric)[condition == "groupA"], get(metric)[condition == "control"])$p.value)
})
wilcoxon_p_df <- data.frame(diversity = diversity_list, p_value = wilcoxon_p)
write.csv(wilcoxon_p_df, file.path(folder, "wilcoxon_p_groupA_vs_control.csv"))

###############################################################################
## Beta-diversity: PERMANOVA and PCoA

# Arcsine sqrt transform for beta diversity
transformed <- asin(sqrt(groupA_vs_control_data[, 1:(ncol(groupA_vs_control_data)-33)]))
permanova_df <- merge(transformed, metadata, by = "row.names", all = FALSE)
rownames(permanova_df) <- permanova_df$Row.names
permanova_df <- permanova_df[, -1]

# PERMANOVA: marginal (Type III-like) effects
set.seed(10)
results_all_condition <- adonis2(
  permanova_df[, 1:(ncol(permanova_df)-32)] ~ age + sex + BMI + condition,
  data = permanova_df,
  permutations = 999,
  strata = permanova_df$household_id,
  by = "margin"
)
capture.output(results_all_condition, file = file.path(folder, "PERMANOVA_results_all_marginal_groupA_vs_control.txt"))

# PERMANOVA: univariate models
set.seed(10)
results_condition_only <- adonis2(permanova_df[, 1:(ncol(permanova_df)-32)] ~ condition, data = permanova_df, permutations = 999, strata = permanova_df$household_id)
capture.output(results_condition_only, file = file.path(folder, "PERMANOVA_results_condition_only_groupA_vs_control.txt"))

set.seed(10)
results_age_only <- adonis2(permanova_df[, 1:(ncol(permanova_df)-32)] ~ age, data = permanova_df, permutations = 999, strata = permanova_df$household_id)
capture.output(results_age_only, file = file.path(folder, "PERMANOVA_results_age_only_groupA_vs_control.txt"))

set.seed(10)
results_BMI_only <- adonis2(permanova_df[, 1:(ncol(permanova_df)-32)] ~ BMI, data = permanova_df, permutations = 999, strata = permanova_df$household_id)
capture.output(results_BMI_only, file = file.path(folder, "PERMANOVA_results_BMI_only_groupA_vs_control.txt"))

set.seed(10)
results_sex_only <- adonis2(permanova_df[, 1:(ncol(permanova_df)-32)] ~ sex, data = permanova_df, permutations = 999, strata = permanova_df$household_id)
capture.output(results_sex_only, file = file.path(folder, "PERMANOVA_results_sex_only_groupA_vs_control.txt"))

# PCoA plot
distance_matrix <- vegdist(permanova_df[, 1:(ncol(permanova_df)-32)], method = "bray")
pcoa_res <- dudi.pco(distance_matrix, scannf = FALSE, nf = 3)

# Extract % variance explained
evals <- eigenvals(pcoa_res)
variance <- evals / sum(evals)
pc1_var <- 100 * signif(variance[1], 2)
pc2_var <- 100 * signif(variance[2], 2)

# Format PCoA coordinates
pc_data <- data.frame(
  pc1 = pcoa_res$li$A1,
  pc2 = pcoa_res$li$A2,
  Group = as.factor(permanova_df$condition)
)

# Compute group centroids
centroids <- aggregate(cbind(pc1, pc2) ~ Group, data = pc_data, FUN = mean)
merged_data <- merge(pc_data, centroids, by = "Group")

# Generate PCoA plot
pdf(file.path(folder, "pcoa_groupA_vs_control.pdf"), width = 6, height = 4.5)
ggplot(merged_data, aes(x = c_pc1, y = c_pc2, color = Group)) +
  scale_colour_manual(values = c("#ff9274", "#55b7e6")) +
  theme(legend.title = element_blank()) +
  labs(x = paste("PCoA1 (", pc1_var, "%)"), y = paste("PCoA2 (", pc2_var, "%)")) +
  geom_segment(aes(x = c_pc1, y = c_pc2, xend = pc1, yend = pc2)) +
  geom_point(aes(x = pc1, y = pc2), size = 3, shape = 20) +
  stat_ellipse(aes(x = pc1, y = pc2, group = Group), level = 0.95) +
  theme_bw() +
  theme(panel.grid.major = element_line(linetype = "dashed", color = "gray70", linewidth = 0.3),
        panel.grid.minor = element_blank())
dev.off()
