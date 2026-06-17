# ==========================================================
# DELIMITATION: Refined Fixed Threshold (Meier et al., 2006)
# ==========================================================

# 1. Directory setup
DATA_DIR <- "C:/GEN/Bosmina/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if(!dir.exists(OUT_DIR))  dir.create(OUT_DIR, recursive = TRUE)
setwd(DATA_DIR)

FASTA_FILE <- "_bosmina_all_GB_MEG_DJT.fasta"
TREE_FILE  <- "_bosmina_all_GB_MEG_DJT.fasta.treefile"
OUT_CSV_PATRI <- file.path(OUT_DIR, "Delimitation_Refined_patri.csv")
OUT_CSV_PAIR  <- file.path(OUT_DIR, "Delimitation_Refined_pair.csv")

THRESHOLD_K80    <- 0.03    # 3% for K80
THRESHOLD_PATRI  <- 0.10    # adjusted for patristic distances

# 2. Packages
pkgs <- c("ape")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting Refined Fixed Threshold...\n")

# 3. Reading data
tree <- read.tree(TREE_FILE)
dna  <- read.dna(FASTA_FILE, format = "fasta")
cat(sprintf("Loaded: %d taxa (tree), %d sequences (FASTA)\n",
            length(tree$tip.label), nrow(dna)))

# 4. Outgroup filtering
OUTGROUP_PATTERN <- "Bosminopsis"
outgroup_in_tree <- tree$tip.label[grepl(OUTGROUP_PATTERN, tree$tip.label, ignore.case = TRUE)]
outgroup_in_dna  <- rownames(dna)[grepl(OUTGROUP_PATTERN, rownames(dna), ignore.case = TRUE)]
if(length(outgroup_in_tree) > 0) tree <- drop.tip(tree, outgroup_in_tree)
if(length(outgroup_in_dna) > 0)  dna <- dna[!rownames(dna) %in% outgroup_in_dna, , drop = FALSE]

# 5. Name synchronization
rownames(dna) <- gsub("\\|", "_", rownames(dna))
common_tips <- intersect(tree$tip.label, rownames(dna))
if(length(common_tips) < length(tree$tip.label)) {
  tree <- drop.tip(tree, setdiff(tree$tip.label, common_tips))
}
dna <- dna[common_tips, , drop = FALSE]
cat("Taxa for analysis:", nrow(dna), "\n")

# =====================================================================
# Refined threshold clustering function (single-linkage + validation)
# =====================================================================
refined_threshold_clustering <- function(dist_mat, threshold) {
  mat <- as.matrix(dist_mat)
  tip_names <- rownames(mat)
  
  # Step 1: Single-linkage clustering at fixed threshold
  # This merges all sequences connected by a chain of d < threshold
  hc <- hclust(as.dist(mat), method = "single")
  clusters <- cutree(hc, h = threshold)
  names(clusters) <- tip_names
  
  # Step 2: Refinement — checking intra-cluster distances
  # If max within cluster > threshold, mark cluster as "problematic"
  # and recursively split it via UPGMA at a sub-threshold
  cl_sizes <- table(clusters)
  multi_clusters <- as.integer(names(cl_sizes[cl_sizes > 1]))
  
  for(cl in multi_clusters) {
    idx <- which(clusters == cl)
    sub_mat <- mat[idx, idx, drop = FALSE]
    
    # Checking if there are pairs within the cluster with distance > threshold
    max_intra <- max(sub_mat[upper.tri(sub_mat)])
    
    if(max_intra > threshold * 1.5) {
      # Problematic cluster: split via UPGMA with stricter threshold
      hc_sub <- hclust(as.dist(sub_mat), method = "average")
      
      # Finding optimal sub-threshold (minimum between threshold/2 and threshold)
      heights <- hc_sub$height
      valid_heights <- heights[heights > threshold/2 & heights < threshold]
      
      if(length(valid_heights) > 0) {
        sub_threshold <- median(valid_heights)
      } else {
        sub_threshold <- threshold / 2
      }
      
      sub_clusters <- cutree(hc_sub, h = sub_threshold)
      
      # Renumbering sub-clusters to avoid overlap with main clusters
      max_cl <- max(clusters, na.rm = TRUE)
      for(j in seq_along(sub_clusters)) {
        if(sub_clusters[j] > 1) {
          clusters[idx[j]] <- max_cl + sub_clusters[j] - 1
        }
      }
    }
  }
  
  # Renumbering clusters sequentially
  unique_cl <- unique(clusters)
  mapping <- setNames(seq_along(unique_cl), unique_cl)
  clusters <- as.integer(mapping[as.character(clusters)])
  names(clusters) <- tip_names
  
  return(clusters)
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

cluster_patri <- refined_threshold_clustering(dist_mat_patri_clean, THRESHOLD_PATRI)
cluster_patri <- cluster_patri[tree$tip.label]

cat(sprintf("Threshold (patristic): %.4f\n", THRESHOLD_PATRI))
cat(sprintf("Refined Threshold (patristic): %d MOTUs\n", length(unique(cluster_patri))))

results_patri <- data.frame(
  Sequence     = names(cluster_patri),
  MOTU_Refined = as.character(cluster_patri),
  stringsAsFactors = FALSE
)
results_patri <- results_patri[order(as.numeric(results_patri$MOTU_Refined), results_patri$Sequence), ]
rownames(results_patri) <- NULL

write.table(results_patri, file = OUT_CSV_PATRI,
            row.names = FALSE, sep = ";", dec = ".",
            quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s\n", basename(OUT_CSV_PATRI)))

# =====================================================================
# VARIANT 2: K80 distances
# =====================================================================
cat("\nVARIANT 2: K80 distances...\n")
dist_mat_pair <- dist.dna(dna, model = "K80", pairwise.deletion = TRUE)
dist_mat_pair_clean <- as.dist(as.matrix(dist_mat_pair))

if(any(is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean))) {
  dist_mat_pair_clean[is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean)] <- 0.0001
}

na_ratio <- 1 - sum(!is.na(as.vector(dist_mat_pair))) / length(as.vector(dist_mat_pair))
if(na_ratio > 0) {
  cat(sprintf("Warning: %.1f%% distances = NA (gaps in alignment)\n", na_ratio * 100))
}

cluster_pair <- refined_threshold_clustering(dist_mat_pair_clean, THRESHOLD_K80)
cluster_pair <- cluster_pair[rownames(dna)]

cat(sprintf("Threshold (K80): %.4f\n", THRESHOLD_K80))
cat(sprintf("Refined Threshold (K80): %d MOTUs\n", length(unique(cluster_pair))))

results_pair <- data.frame(
  Sequence     = names(cluster_pair),
  MOTU_Refined = as.character(cluster_pair),
  stringsAsFactors = FALSE
)
results_pair <- results_pair[order(as.numeric(results_pair$MOTU_Refined), results_pair$Sequence), ]
rownames(results_pair) <- NULL

write.table(results_pair, file = OUT_CSV_PAIR,
            row.names = FALSE, sep = ";", dec = ".",
            quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s\n", basename(OUT_CSV_PAIR)))

# =====================================================================
# Summary with cluster sizes
# =====================================================================
cl_sizes_patri <- table(cluster_patri)
cl_sizes_pair  <- table(cluster_pair)

cat("\nSummary:\n")
cat("   Patristic:\n")
cat(sprintf("      - Threshold:          %.4f\n", THRESHOLD_PATRI))
cat(sprintf("      - Total MOTUs:        %d\n", length(cl_sizes_patri)))
cat(sprintf("      - Cluster sizes:      from %d to %d (median %.0f)\n",
            min(cl_sizes_patri), max(cl_sizes_patri), median(cl_sizes_patri)))
cat(sprintf("      - Singletons:         %d\n", sum(cl_sizes_patri == 1)))

cat("   K80:\n")
cat(sprintf("      - Threshold:          %.4f\n", THRESHOLD_K80))
cat(sprintf("      - Total MOTUs:        %d\n", length(cl_sizes_pair)))
cat(sprintf("      - Cluster sizes:      from %d to %d (median %.0f)\n",
            min(cl_sizes_pair), max(cl_sizes_pair), median(cl_sizes_pair)))
cat(sprintf("      - Singletons:         %d\n", sum(cl_sizes_pair == 1)))

cat(sprintf("\nRefined Fixed Threshold completed. Results in:\n   %s\n", OUT_DIR))