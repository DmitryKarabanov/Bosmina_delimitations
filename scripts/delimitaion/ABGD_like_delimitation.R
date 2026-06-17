# ==========================================================
# DELIMITATION: ABGD-like (recursive partitioning, local implementation)
# ==========================================================

# 1. Directory setup
DATA_DIR <- "C:/GEN/Bosmina/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if(!dir.exists(OUT_DIR))  dir.create(OUT_DIR, recursive = TRUE)
setwd(DATA_DIR)

FASTA_FILE <- "_bosmina_all_GB_MEG_DJT.fasta"
TREE_FILE  <- "_bosmina_all_GB_MEG_DJT.fasta.treefile"

OUT_CSV_PATRI_INIT  <- file.path(OUT_DIR, "Delimitation_ABGD_init_patri.csv")
OUT_CSV_PATRI_RECUR <- file.path(OUT_DIR, "Delimitation_ABGD_recur_patri.csv")
OUT_CSV_PAIR_INIT   <- file.path(OUT_DIR, "Delimitation_ABGD_init_pair.csv")
OUT_CSV_PAIR_RECUR  <- file.path(OUT_DIR, "Delimitation_ABGD_recur_pair.csv")

# 2. Packages
pkgs <- c("ape")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting ABGD-like analysis...\n")

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
cat("Taxa for analysis:", nrow(dna), "\n")

# =====================================================================
# ABGD-like algorithm core
# =====================================================================

# Barcode gap search in distance distribution (local density minimum)
find_barcode_gap <- function(dist_vec) {
  if(length(dist_vec) < 3) return(NA)
  den <- density(dist_vec, from = 0, n = 256, adjust = 1.5)
  y <- den$y
  x <- den$x
  minima <- c()
  
  for(i in 2:(length(y) - 1)) {
    if(y[i] < y[i - 1] && y[i] < y[i + 1]) {
      # Minimum must lie in a reasonable range (not in the tails)
      if(x[i] > 0.005 && x[i] < max(dist_vec) * 0.5) {
        minima <- c(minima, x[i])
      }
    }
  }
  if(length(minima) == 0) return(NA)
  return(minima[1])  # first suitable minimum
}

# Recursive partitioning (reproduces the logic of 'recursive partition' of the original ABGD)
recursive_partition <- function(dist_mat, indices, depth = 0, max_depth = 5) {
  if(length(indices) < 3 || depth >= max_depth) {
    return(list(indices))
  }
  
  sub_mat <- as.matrix(dist_mat)[indices, indices, drop = FALSE]
  dist_vec <- sub_mat[lower.tri(sub_mat)]
  dist_vec <- dist_vec[!is.na(dist_vec) & dist_vec > 0]
  
  if(length(dist_vec) < 3) return(list(indices))
  
  gap <- find_barcode_gap(dist_vec)
  if(is.na(gap)) return(list(indices))
  
  hc <- hclust(as.dist(sub_mat), method = "average")
  clusters <- cutree(hc, h = gap)
  unique_cl <- unique(clusters)
  
  if(length(unique_cl) <= 1) return(list(indices))
  
  # Recursively partitioning each subcluster
  result <- list()
  for(cl in unique_cl) {
    sub_idx <- indices[which(clusters == cl)]
    if(length(sub_idx) >= 3 && depth < max_depth) {
      sub_parts <- recursive_partition(dist_mat, sub_idx, depth + 1, max_depth)
      result <- c(result, sub_parts)
    } else {
      result <- c(result, list(sub_idx))
    }
  }
  return(result)
}

# Main ABGD-like execution function
run_abgd_like <- function(dist_mat, tip_labels, partition_type = "initial") {
  n <- length(tip_labels)
  
  # Initial partition: single pass barcode gap search
  dist_vec <- as.vector(as.matrix(dist_mat)[lower.tri(dist_mat)])
  dist_vec <- dist_vec[!is.na(dist_vec)]
  initial_gap <- find_barcode_gap(dist_vec)
  
  if(is.na(initial_gap)) {
    cat("Barcode gap not found, using distance median as fallback\n")
    initial_gap <- median(dist_vec)
  }
  
  hc <- hclust(dist_mat, method = "average")
  initial_cl <- cutree(hc, h = initial_gap)
  
  if(partition_type == "initial") {
    return(initial_cl)
  }
  
  # Recursive partition: recursive partitioning
  parts <- recursive_partition(dist_mat, 1:n, depth = 0, max_depth = 5)
  recursive_cl <- rep(NA, n)
  names(recursive_cl) <- tip_labels
  
  for(i in seq_along(parts)) {
    recursive_cl[parts[[i]]] <- i
  }
  return(recursive_cl)
}

# Result saving function
save_results <- function(cluster_vec, out_file) {
  results <- data.frame(
    Sequence  = names(cluster_vec),
    MOTU_ABGD = as.character(cluster_vec),
    stringsAsFactors = FALSE
  )
  results <- results[order(as.numeric(results$MOTU_ABGD), results$Sequence), ]
  rownames(results) <- NULL
  
  write.table(results, file = out_file, row.names = FALSE, sep = ";", dec = ".",
              quote = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("Saved: %s (%d MOTUs)\n", basename(out_file), length(unique(cluster_vec))))
}

# =====================================================================
# VARIANT 1: Patristic distances
# =====================================================================
cat("\nVARIANT 1: Patristic distances...\n")
dist_mat_patri <- cophenetic(tree)
dist_mat_patri_clean <- as.dist(as.matrix(dist_mat_patri))

if(any(is.na(dist_mat_patri_clean) | is.nan(dist_mat_patri_clean))) {
  dist_mat_patri_clean[is.na(dist_mat_patri_clean) | is.nan(dist_mat_patri_clean)] <- 0.0001
}

cluster_patri_init <- run_abgd_like(dist_mat_patri_clean, tree$tip.label, "initial")
cluster_patri_init <- cluster_patri_init[tree$tip.label]
save_results(cluster_patri_init, OUT_CSV_PATRI_INIT)

cluster_patri_recur <- run_abgd_like(dist_mat_patri_clean, tree$tip.label, "recursive")
cluster_patri_recur <- cluster_patri_recur[tree$tip.label]
save_results(cluster_patri_recur, OUT_CSV_PATRI_RECUR)

# =====================================================================
# VARIANT 2: Pairwise K80 distances
# =====================================================================
cat("\nVARIANT 2: Pairwise K80 distances...\n")
dist_mat_pair <- dist.dna(dna, model = "K80", pairwise.deletion = TRUE)
dist_mat_pair_clean <- as.dist(as.matrix(dist_mat_pair))

if(any(is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean))) {
  dist_mat_pair_clean[is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean)] <- 0.0001
}

na_ratio <- 1 - sum(!is.na(as.vector(dist_mat_pair))) / length(as.vector(dist_mat_pair))
if(na_ratio > 0) {
  cat(sprintf("Warning: %.1f%% distances = NA (gaps in alignment)\n", na_ratio * 100))
}

cluster_pair_init <- run_abgd_like(dist_mat_pair_clean, rownames(dna), "initial")
cluster_pair_init <- cluster_pair_init[rownames(dna)]
save_results(cluster_pair_init, OUT_CSV_PAIR_INIT)

cluster_pair_recur <- run_abgd_like(dist_mat_pair_clean, rownames(dna), "recursive")
cluster_pair_recur <- cluster_pair_recur[rownames(dna)]
save_results(cluster_pair_recur, OUT_CSV_PAIR_RECUR)

# =====================================================================
# Agreement assessment between initial and recursive partitions
# =====================================================================
cat("\nEvaluating agreement between initial and recursive partitions...\n")
calc_agreement <- function(cl1, cl2) {
  valid <- !(is.na(cl1) | is.na(cl2))
  if(sum(valid) < 2) return(NA)
  cl1 <- cl1[valid]; cl2 <- cl2[valid]
  s1 <- outer(cl1, cl1, "=="); s2 <- outer(cl2, cl2, "==")
  agree <- (s1 & s2) | (!s1 & !s2)
  n <- length(cl1)
  return(sum(agree[upper.tri(agree)]) / (n * (n - 1) / 2) * 100)
}

agr_patri <- calc_agreement(cluster_patri_init[tree$tip.label], cluster_patri_recur[tree$tip.label])
agr_pair  <- calc_agreement(cluster_pair_init[rownames(dna)], cluster_pair_recur[rownames(dna)])

# =====================================================================
# Final summary
# =====================================================================
cat("\nSummary:\n")
cat("   Patristic:\n")
cat(sprintf("      - Initial partition:   %d MOTUs\n", length(unique(cluster_patri_init))))
cat(sprintf("      - Recursive partition: %d MOTUs\n", length(unique(cluster_patri_recur))))
cat(sprintf("      - Agreement:           %.1f%%\n", agr_patri))

cat("   K80:\n")
cat(sprintf("      - Initial partition:   %d MOTUs\n", length(unique(cluster_pair_init))))
cat(sprintf("      - Recursive partition: %d MOTUs\n", length(unique(cluster_pair_recur))))
cat(sprintf("      - Agreement:           %.1f%%\n", agr_pair))

cat(sprintf("\nABGD-like analysis completed. Results in:\n   %s\n", OUT_DIR))