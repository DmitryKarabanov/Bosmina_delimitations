# ==========================================================
# SPECIES DELIMITATION: locMin (two distance variants)
# ==========================================================

# 1. Directory setup
DATA_DIR <- "C:/GEN/Bosmina/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if(!dir.exists(OUT_DIR))  dir.create(OUT_DIR, recursive = TRUE)
setwd(DATA_DIR)

FASTA_FILE <- "_bosmina_all_GB_MEG_DJT.fasta"
TREE_FILE  <- "_bosmina_all_GB_MEG_DJT.fasta.treefile"
OUT_CSV_PATRI <- file.path(OUT_DIR, "Delimitation_locMin_patri.csv")
OUT_CSV_PAIR  <- file.path(OUT_DIR, "Delimitation_locMin_pair.csv")
OUT_PNG       <- file.path(OUT_DIR, "Distances_histogram_locMin_both.png")

# 2. Packages
pkgs <- c("ape", "ggplot2")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting locMin analysis (2 variants)...\n")

# 3. Reading data
tree <- read.tree(TREE_FILE)
dna  <- read.dna(FASTA_FILE, format = "fasta")
cat(sprintf("Loaded: %d taxa (tree), %d sequences (FASTA)\n",
            length(tree$tip.label), nrow(dna)))

# 4. Outgroup filtering (Bosminopsis)
OUTGROUP_PATTERN <- "Bosminopsis"
outgroup_in_tree <- tree$tip.label[grepl(OUTGROUP_PATTERN, tree$tip.label, ignore.case = TRUE)]
outgroup_in_dna  <- rownames(dna)[grepl(OUTGROUP_PATTERN, rownames(dna), ignore.case = TRUE)]
if(length(outgroup_in_tree) > 0) tree <- drop.tip(tree, outgroup_in_tree)
if(length(outgroup_in_dna) > 0)  dna <- dna[!rownames(dna) %in% outgroup_in_dna, , drop = FALSE]
cat(sprintf("Removed outgroups: %d taxa\n", length(unique(c(outgroup_in_tree, outgroup_in_dna)))))

# 5. Name synchronization (tree and FASTA must have identical names)
rownames(dna) <- gsub("\\|", "_", rownames(dna))
common_tips <- intersect(tree$tip.label, rownames(dna))
if(length(common_tips) < length(tree$tip.label)) {
  tree <- drop.tip(tree, setdiff(tree$tip.label, common_tips))
}
dna <- dna[common_tips, , drop = FALSE]
cat("Taxa for analysis:", nrow(dna), "\n")

# =====================================================================
# VARIANT 1: Patristic distances (from IQ-TREE)
# =====================================================================
cat("\nVARIANT 1: Patristic distances...\n")
dist_mat_patri <- cophenetic(tree)
dist_vec_patri <- as.vector(dist_mat_patri[lower.tri(dist_mat_patri)])

localMinima2 <- function(distobj) {
  den <- density(distobj, from = 0, n = 512, adjust = 1.2)
  a <- rep(NA, length(den$y) - 2)
  for (i in 2:(length(den$y) - 1)) {
    a[i - 1] <- den$y[i - 1] > den$y[i] & den$y[i + 1] > den$y[i]
  }
  den$localMinima <- den$x[which(a)]
  invisible(den)
}

lmin_patri <- localMinima2(dist_vec_patri)
threshold_patri <- if(length(lmin_patri$localMinima) > 0) lmin_patri$localMinima[1] else 0.03
cat(sprintf("locMin threshold (patristic): %.4f\n", threshold_patri))

dist_mat_patri_clean <- as.dist(as.matrix(dist_mat_patri))
if(any(is.na(dist_mat_patri_clean) | is.nan(dist_mat_patri_clean))) {
  dist_mat_patri_clean[is.na(dist_mat_patri_clean) | is.nan(dist_mat_patri_clean)] <- 0.0001
}

hc_patri <- hclust(dist_mat_patri_clean, method = "average")
cluster_patri <- cutree(hc_patri, h = threshold_patri)
cluster_patri <- cluster_patri[tree$tip.label]

results_patri <- data.frame(
  Sequence    = names(cluster_patri),
  MOTU_locMin = as.character(cluster_patri),
  stringsAsFactors = FALSE
)
results_patri <- results_patri[order(as.numeric(results_patri$MOTU_locMin), results_patri$Sequence), ]
rownames(results_patri) <- NULL

write.table(results_patri, file = OUT_CSV_PATRI,
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s (%d rows, %d MOTUs)\n",
            OUT_CSV_PATRI, nrow(results_patri), length(unique(cluster_patri))))

# =====================================================================
# VARIANT 2: Pairwise K80 distances (from FASTA)
# =====================================================================
cat("\nVARIANT 2: Pairwise K80 distances...\n")
dist_mat_pair <- dist.dna(dna, model = "K80", pairwise.deletion = TRUE)
dist_vec_pair <- as.vector(dist_mat_pair[lower.tri(dist_mat_pair)])
dist_vec_pair_clean <- dist_vec_pair[!is.na(dist_vec_pair)]
na_ratio <- 1 - length(dist_vec_pair_clean) / length(dist_vec_pair)

if(na_ratio > 0) {
  cat(sprintf("Warning: %.1f%% distances = NA (gaps in alignment)\n", na_ratio * 100))
}

lmin_pair <- localMinima2(dist_vec_pair_clean)
threshold_pair <- if(length(lmin_pair$localMinima) > 0) lmin_pair$localMinima[1] else 0.03
cat(sprintf("locMin threshold (K80): %.4f\n", threshold_pair))

dist_mat_pair_clean <- as.dist(as.matrix(dist_mat_pair))
if(any(is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean))) {
  dist_mat_pair_clean[is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean)] <- 0.0001
}

hc_pair <- hclust(dist_mat_pair_clean, method = "average")
cluster_pair <- cutree(hc_pair, h = threshold_pair)
cluster_pair <- cluster_pair[rownames(dna)]

results_pair <- data.frame(
  Sequence    = names(cluster_pair),
  MOTU_locMin = as.character(cluster_pair),
  stringsAsFactors = FALSE
)
results_pair <- results_pair[order(as.numeric(results_pair$MOTU_locMin), results_pair$Sequence), ]
rownames(results_pair) <- NULL

write.table(results_pair, file = OUT_CSV_PAIR,
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s (%d rows, %d MOTUs)\n",
            OUT_CSV_PAIR, nrow(results_pair), length(unique(cluster_pair))))

# =====================================================================
# Comparative visualization (two histograms with thresholds)
# =====================================================================
df_patri <- data.frame(dist = dist_vec_patri, type = "Patristic (IQ-TREE)")
df_pair  <- data.frame(dist = dist_vec_pair_clean, type = "K80 (FASTA)")
df_both  <- rbind(df_patri, df_pair)

p_hist <- ggplot(df_both, aes(x = dist, fill = type)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity", color = "white") +
  geom_vline(data = data.frame(type = c("Patristic (IQ-TREE)", "K80 (FASTA)"),
                               thr = c(threshold_patri, threshold_pair)),
             aes(xintercept = thr, color = type), linewidth = 1.2, linetype = "dashed") +
  scale_fill_manual(values = c("Patristic (IQ-TREE)" = "steelblue", "K80 (FASTA)" = "coral")) +
  scale_color_manual(values = c("Patristic (IQ-TREE)" = "#DC143C", "K80 (FASTA)" = "#228B22")) +
  facet_wrap(~type, ncol = 1, scales = "free_y") +
  labs(title = "locMin: comparison of two distance metrics",
       subtitle = paste("Patristic: threshold =", round(threshold_patri, 4), "|",
                        "K80: threshold =", round(threshold_pair, 4)),
       x = "Genetic distance", y = "Frequency") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

ggsave(filename = OUT_PNG, plot = p_hist, width = 10, height = 8, dpi = 300)

# =====================================================================
# Agreement assessment between the two variants
# =====================================================================
cat("\nEvaluating agreement between the two locMin variants...\n")

# Synchronizing order
cluster_patri_aligned <- cluster_patri[tree$tip.label]
cluster_pair_aligned  <- cluster_pair[tree$tip.label]

calc_agreement <- function(cl1, cl2) {
  valid <- !(is.na(cl1) | is.na(cl2))
  if(sum(valid) < 2) return(NA)
  cl1 <- cl1[valid]; cl2 <- cl2[valid]
  s1 <- outer(cl1, cl1, "=="); s2 <- outer(cl2, cl2, "==")
  agree <- (s1 & s2) | (!s1 & !s2)
  n <- length(cl1)
  return(sum(agree[upper.tri(agree)]) / (n*(n-1)/2) * 100)
}

agreement <- calc_agreement(cluster_patri_aligned, cluster_pair_aligned)

# =====================================================================
# Final summary
# =====================================================================
cl_sizes_patri <- table(cluster_patri)
cl_sizes_pair  <- table(cluster_pair)

cat("\nSummary:\n")
cat("   Variant 1 (Patristic):\n")
cat(sprintf("      - Threshold:          %.4f\n", threshold_patri))
cat(sprintf("      - Total MOTUs:        %d\n", length(cl_sizes_patri)))
cat(sprintf("      - Cluster sizes:      from %d to %d (median %.0f)\n",
            min(cl_sizes_patri), max(cl_sizes_patri), median(cl_sizes_patri)))
cat(sprintf("      - Singletons:         %d\n", sum(cl_sizes_patri == 1)))

cat("   Variant 2 (K80):\n")
cat(sprintf("      - Threshold:          %.4f\n", threshold_pair))
cat(sprintf("      - Total MOTUs:        %d\n", length(cl_sizes_pair)))
cat(sprintf("      - Cluster sizes:      from %d to %d (median %.0f)\n",
            min(cl_sizes_pair), max(cl_sizes_pair), median(cl_sizes_pair)))
cat(sprintf("      - Singletons:         %d\n", sum(cl_sizes_pair == 1)))

cat(sprintf("\nMethod agreement: %.1f%% of sequence pairs grouped identically\n", agreement))

if(agreement >= 80) {
  cat("Conclusion: both variants yield consistent results. Delimitation is robust.\n")
} else if(agreement >= 60) {
  cat("Conclusion: moderate agreement. K80 variant splits more strongly (expected).\n")
} else {
  cat("Conclusion: strong discrepancy. Check alignment quality and sequence lengths.\n")
}

cat(sprintf("\nAnalysis completed. Results in:\n   %s\n", OUT_DIR))
cat("   Delimitation_locMin_patri.csv (patristic)\n")
cat("   Delimitation_locMin_pair.csv (K80)\n")
cat("   Distances_histogram_locMin_both.png\n")