# ==========================================================
# DELIMITATION: Ward's Method (two distance variants)
# ==========================================================

# 1. Directory setup
DATA_DIR <- "C:/GEN/Bosmina/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)
if(!dir.exists(OUT_DIR))  dir.create(OUT_DIR, recursive = TRUE)
setwd(DATA_DIR)

FASTA_FILE <- "_bosmina_all_GB_MEG_DJT.fasta"
TREE_FILE  <- "_bosmina_all_GB_MEG_DJT.fasta.treefile"
OUT_CSV_PATRI <- file.path(OUT_DIR, "Delimitation_Ward_patri.csv")
OUT_CSV_PAIR  <- file.path(OUT_DIR, "Delimitation_Ward_pair.csv")

# 2. Packages
pkgs <- c("ape")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting Ward's Method...\n")

# 3. Reading data
tree <- read.tree(TREE_FILE)
dna  <- read.dna(FASTA_FILE, format = "fasta")

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
# VARIANT 1: Patristic distances
# =====================================================================
cat("\nVARIANT 1: Patristic distances...\n")
dist_mat_patri <- cophenetic(tree)
dist_vec_patri <- as.vector(dist_mat_patri[lower.tri(dist_mat_patri)])
dist_mat_patri_clean <- as.dist(as.matrix(dist_mat_patri))

if(any(is.na(dist_mat_patri_clean) | is.nan(dist_mat_patri_clean))) {
  dist_mat_patri_clean[is.na(dist_mat_patri_clean) | is.nan(dist_mat_patri_clean)] <- 0.0001
}

# Ward's method clustering
hc_patri <- hclust(dist_mat_patri_clean, method = "ward.D2")

# locMin threshold search
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
cluster_patri <- cutree(hc_patri, h = threshold_patri)
cluster_patri <- cluster_patri[tree$tip.label]

cat(sprintf("Ward threshold (patristic): %.4f\n", threshold_patri))
cat(sprintf("Ward (patristic): %d MOTUs\n", length(unique(cluster_patri))))

results_patri <- data.frame(
  Sequence    = names(cluster_patri),
  MOTU_Ward   = as.character(cluster_patri),
  stringsAsFactors = FALSE
)
results_patri <- results_patri[order(as.numeric(results_patri$MOTU_Ward), results_patri$Sequence), ]
rownames(results_patri) <- NULL

write.table(results_patri, file = OUT_CSV_PATRI,
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s\n", OUT_CSV_PATRI))

# =====================================================================
# VARIANT 2: K80 distances
# =====================================================================
cat("\nVARIANT 2: K80 distances...\n")
dist_mat_pair <- dist.dna(dna, model = "K80", pairwise.deletion = TRUE)
dist_vec_pair <- as.vector(dist_mat_pair[lower.tri(dist_mat_pair)])
dist_vec_pair_clean <- dist_vec_pair[!is.na(dist_vec_pair)]
dist_mat_pair_clean <- as.dist(as.matrix(dist_mat_pair))

if(any(is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean))) {
  dist_mat_pair_clean[is.na(dist_mat_pair_clean) | is.nan(dist_mat_pair_clean)] <- 0.0001
}

hc_pair <- hclust(dist_mat_pair_clean, method = "ward.D2")
lmin_pair <- localMinima2(dist_vec_pair_clean)
threshold_pair <- if(length(lmin_pair$localMinima) > 0) lmin_pair$localMinima[1] else 0.03
cluster_pair <- cutree(hc_pair, h = threshold_pair)
cluster_pair <- cluster_pair[rownames(dna)]

cat(sprintf("Ward threshold (K80): %.4f\n", threshold_pair))
cat(sprintf("Ward (K80): %d MOTUs\n", length(unique(cluster_pair))))

results_pair <- data.frame(
  Sequence    = names(cluster_pair),
  MOTU_Ward   = as.character(cluster_pair),
  stringsAsFactors = FALSE
)
results_pair <- results_pair[order(as.numeric(results_pair$MOTU_Ward), results_pair$Sequence), ]
rownames(results_pair) <- NULL

write.table(results_pair, file = OUT_CSV_PAIR,
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s\n", OUT_CSV_PAIR))

cat(sprintf("\nWard's Method completed. Results in:\n   %s\n", OUT_DIR))