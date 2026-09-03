#!/usr/bin/env Rscript
# Host-microbe interaction analysis for diurnal saliva pilot
#
# Required packages:
#   install.packages(c("vegan", "lme4", "ggplot2", "dplyr", "tidyr", "readr", "tibble"))

library(vegan)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: diurnal_analysis.R <feature_tables_dir> <output_dir> [wms_species_tsv]")
}
feature_dir <- args[1]
out_dir <- args[2]
wms_path <- if (length(args) >= 3) args[3] else NULL

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
micro <- read_tsv(file.path(feature_dir, "microbiome_clr.tsv"), show_col_types = FALSE) %>%
  column_to_rownames("...1")
host <- read_tsv(file.path(feature_dir, "host_features.tsv"), show_col_types = FALSE) %>%
  column_to_rownames("...1")
meta <- read_tsv(file.path(feature_dir, "sample_metadata.tsv"), show_col_types = FALSE) %>%
  column_to_rownames("...1")

# Ensure same order
common <- intersect(rownames(micro), rownames(meta))
micro <- micro[common, , drop = FALSE]
host <- host[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

# Replace missing host dosages (-1) with mean
host[host == -1] <- NA
host <- as.data.frame(lapply(host, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}))

# Host-genotype PCs for candidate SNPs
snp_cols <- setdiff(colnames(host), "host_fraction")
if (length(snp_cols) > 1) {
  host_pca <- prcomp(host[, snp_cols, drop = FALSE], scale. = TRUE)
  host$host_geno_PC1 <- host_pca$x[, 1]
  host$host_geno_PC2 <- host_pca$x[, 2]
} else {
  host$host_geno_PC1 <- 0
  host$host_geno_PC2 <- 0
}

meta$host_fraction <- host$host_fraction
meta$host_geno_PC1 <- host$host_geno_PC1
meta$host_geno_PC2 <- host$host_geno_PC2

# Beta diversity
# Use raw counts if available; here we use CLR as input for Bray-Curtis approximation
# Better: re-load raw species_counts and proportion-normalize.
raw_counts <- read_tsv(file.path(feature_dir, "microbiome_clr.tsv"), show_col_types = FALSE) %>%
  column_to_rownames("...1")
prop <- raw_counts / rowSums(raw_counts)
prop <- prop[, colSums(prop > 0) >= 3, drop = FALSE]
prop <- prop[rowSums(prop) > 0, , drop = FALSE]
meta <- meta[rownames(prop), , drop = FALSE]
host <- host[rownames(prop), , drop = FALSE]

dist_bc <- vegdist(prop, method = "bray")

# PERMANOVA
write("\n=== PERMANOVA: distance ~ subject + time_point + host_fraction ===", stdout())
perm1 <- adonis2(dist_bc ~ subject_id + time_point + host_fraction, data = meta, by = "margin")
print(perm1)
write_tsv(as.data.frame(perm1), file.path(out_dir, "permanova_subject_time_hostfrac.tsv"))

write("\n=== PERMANOVA: distance ~ subject + time_point + host_geno_PC1 + host_geno_PC2 ===", stdout())
perm2 <- adonis2(dist_bc ~ subject_id + time_point + host_geno_PC1 + host_geno_PC2,
                 data = meta, by = "margin")
print(perm2)
write_tsv(as.data.frame(perm2), file.path(out_dir, "permanova_subject_time_hostgeno.tsv"))

# PCoA
pcoa <- cmdscale(dist_bc, k = 2, eig = TRUE)
meta$PCo1 <- pcoa$points[, 1]
meta$PCo2 <- pcoa$points[, 2]

p <- ggplot(meta, aes(x = PCo1, y = PCo2, color = subject_id, shape = time_point)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCoA of 2bRAD-M species (Bray-Curtis)")
ggsave(file.path(out_dir, "pcoa_subject_time.pdf"), p, width = 7, height = 5)

p2 <- ggplot(meta, aes(x = PCo1, y = PCo2, color = host_fraction)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "PCoA colored by host fraction")
ggsave(file.path(out_dir, "pcoa_host_fraction.pdf"), p2, width = 6, height = 5)

# Per-taxon mixed models
# Use CLR abundance; keep taxa present in >= 25% of samples
clr <- micro[rownames(prop), , drop = FALSE]
taxa_to_test <- names(which(colSums(clr != 0) >= nrow(clr) * 0.25))

mixed_results <- lapply(taxa_to_test, function(taxon) {
  df <- meta %>% mutate(abn = clr[, taxon])
  m <- lmer(abn ~ time_point + host_fraction + (1 | subject_id), data = df)
  s <- summary(m)$coefficients
  data.frame(
    taxon = taxon,
    term = rownames(s),
    estimate = s[, "Estimate"],
    std_error = s[, "Std. Error"],
    p_value = s[, "Pr(>|t|)"]
  )
})
mixed_df <- bind_rows(mixed_results)
mixed_df$p_adj <- p.adjust(mixed_df$p_value, method = "BH")
write_tsv(mixed_df, file.path(out_dir, "mixed_model_taxa.tsv"))

# Optional Procrustes with WMS
if (!is.null(wms_path) && file.exists(wms_path)) {
  wms <- read_tsv(wms_path, show_col_types = FALSE) %>%
    column_to_rownames("...1")
  common_taxa <- intersect(colnames(prop), colnames(wms))
  wms <- wms[rownames(prop), common_taxa, drop = FALSE]
  wms_prop <- wms / rowSums(wms)
  dist_wms <- vegdist(wms_prop, method = "bray")
  pcoa_wms <- cmdscale(dist_wms, k = 2)
  proc <- procrustes(pcoa$points, pcoa_wms)
  write("\n=== Procrustes 2bRAD-M vs WMS ===", stdout())
  print(summary(proc))
  pdf(file.path(out_dir, "procrustes_2bradm_wms.pdf"), width = 6, height = 5)
  plot(proc)
  dev.off()
}

write("\nAnalysis complete.", stdout())
