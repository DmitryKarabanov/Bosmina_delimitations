# ==========================================================
# DELIMITATION: GMYC (single + multi threshold)
# ==========================================================

DATA_DIR <- "C:/GEN/Bosmina/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

TREE_FILE <- file.path(DATA_DIR, "_bosmina_all_GB_MEG_DJT.fasta.treefile")
OUT_CSV_SINGLE <- file.path(OUT_DIR, "Delimitation_GMYC_single.csv")
OUT_CSV_MULTI  <- file.path(OUT_DIR, "Delimitation_GMYC_multi.csv")

# 1. Packages
pkgs <- c("ape", "splits")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting GMYC analysis...\n")

tree <- read.tree(TREE_FILE)
cat(sprintf("Loaded: %d taxa in the tree\n", length(tree$tip.label)))

# 2. Outgroup filtering
OUTGROUP_PATTERN <- "Bosminopsis"
outgroup_tips <- tree$tip.label[grepl(OUTGROUP_PATTERN, tree$tip.label, ignore.case = TRUE)]
if(length(outgroup_tips) > 0) {
  tree.noroot <- drop.tip(tree, outgroup_tips)
  cat(sprintf("Removed outgroups: %d taxa\n", length(outgroup_tips)))
} else {
  tree.noroot <- tree
}
cat("Taxa for analysis:", length(tree.noroot$tip.label), "\n")

# 3. Ultrametrization
cat("\nUltrametrizing tree via ape::chronos()...\n")
tree.noroot$edge.length[tree.noroot$edge.length < 1e-8] <- 1e-8
tree.ultra <- tryCatch(
  ape::chronos(tree.noroot, model = "correlated", control = list(maxit = 1000)),
  error = function(e) {
    cat("chronos did not converge. Using manual tip alignment...\n")
    n_tips <- length(tree.noroot$tip.label)
    root_idx <- n_tips + 1
    dists <- ape::dist.nodes(tree.noroot)[root_idx, 1:n_tips]
    max_d <- max(dists)
    tree_tmp <- tree.noroot
    tip_edges <- which(tree_tmp$edge[, 2] <= n_tips)
    tree_tmp$edge.length[tip_edges] <- tree_tmp$edge.length[tip_edges] + (max_d - dists[tree_tmp$edge[tip_edges, 2]])
    return(tree_tmp)
  }
)
if(ape::is.ultrametric(tree.ultra)) cat("Tree successfully ultrametrized\n")

# =====================================================================
# VARIANT 1: GMYC single-threshold
# =====================================================================
cat("\nRunning GMYC (single-threshold)...\n")
gmyc_single <- tryCatch(splits::gmyc(tree.ultra, method = "s"),
                        error = function(e) { cat(sprintf("Error: %s\n", e$message)); NULL })

if(!is.null(gmyc_single)) {
  cat("\nGMYC Single-Threshold Summary:\n")
  print(summary(gmyc_single))
  
  spp_single <- spec.list(gmyc_single)
  gmyc_df_single <- data.frame(
    Sequence  = as.character(spp_single$sample_name),
    MOTU_GMYC = as.integer(spp_single$GMYC_spec),
    stringsAsFactors = FALSE
  )
  gmyc_df_single <- gmyc_df_single[order(gmyc_df_single$MOTU_GMYC, gmyc_df_single$Sequence), ]
  rownames(gmyc_df_single) <- NULL
  
  write.table(gmyc_df_single, file = OUT_CSV_SINGLE, row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
  
  n_motu_single <- length(unique(gmyc_df_single$MOTU_GMYC))
  cl_sizes_single <- table(gmyc_df_single$MOTU_GMYC)
  
  cat(sprintf("\nGMYC single: %d MOTUs\n", n_motu_single))
  cat(sprintf("Table saved: %s\n", basename(OUT_CSV_SINGLE)))
} else {
  n_motu_single <- NA
  cl_sizes_single <- NULL
}

# =====================================================================
# VARIANT 2: GMYC multi-threshold
# =====================================================================
cat("\nRunning GMYC (multi-threshold)...\n")
gmyc_multi <- tryCatch(splits::gmyc(tree.ultra, method = "m"),
                       error = function(e) { cat(sprintf("Error: %s\n", e$message)); NULL })

if(!is.null(gmyc_multi)) {
  cat("\nGMYC Multi-Threshold Summary:\n")
  print(summary(gmyc_multi))
  
  spp_multi <- spec.list(gmyc_multi)
  gmyc_df_multi <- data.frame(
    Sequence  = as.character(spp_multi$sample_name),
    MOTU_GMYC = as.integer(spp_multi$GMYC_spec),
    stringsAsFactors = FALSE
  )
  gmyc_df_multi <- gmyc_df_multi[order(gmyc_df_multi$MOTU_GMYC, gmyc_df_multi$Sequence), ]
  rownames(gmyc_df_multi) <- NULL
  
  write.table(gmyc_df_multi, file = OUT_CSV_MULTI, row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
  
  n_motu_multi <- length(unique(gmyc_df_multi$MOTU_GMYC))
  cl_sizes_multi <- table(gmyc_df_multi$MOTU_GMYC)
  
  cat(sprintf("\nGMYC multi: %d MOTUs\n", n_motu_multi))
  cat(sprintf("Table saved: %s\n", basename(OUT_CSV_MULTI)))
} else {
  n_motu_multi <- NA
  cl_sizes_multi <- NULL
}

# =====================================================================
# Comparison single vs multi
# =====================================================================
if(!is.null(gmyc_single) && !is.null(gmyc_multi)) {
  calc_agreement <- function(cl1, cl2) {
    valid <- !(is.na(cl1) | is.na(cl2))
    if(sum(valid) < 2) return(NA)
    cl1 <- cl1[valid]; cl2 <- cl2[valid]
    s1 <- outer(cl1, cl1, "=="); s2 <- outer(cl2, cl2, "==")
    agree <- (s1 & s2) | (!s1 & !s2)
    n <- length(cl1)
    return(sum(agree[upper.tri(agree)]) / (n*(n-1)/2) * 100)
  }
  agreement <- calc_agreement(gmyc_df_single$MOTU_GMYC, gmyc_df_multi$MOTU_GMYC)
  cat(sprintf("\nAgreement single vs multi: %.1f%%\n", agreement))
}

# =====================================================================
# Final summary
# =====================================================================
cat("\nFinal summary:\n")
if(!is.null(gmyc_single)) {
  cat("\n   Single-threshold:\n")
  cat(sprintf("      - Total MOTUs:        %d\n", n_motu_single))
  cat(sprintf("      - Cluster sizes:      from %d to %d (median %.0f)\n", min(cl_sizes_single), max(cl_sizes_single), median(cl_sizes_single)))
  cat(sprintf("      - Singletons:         %d\n", sum(cl_sizes_single == 1)))
}
if(!is.null(gmyc_multi)) {
  cat("\n   Multi-threshold:\n")
  cat(sprintf("      - Total MOTUs:        %d\n", n_motu_multi))
  cat(sprintf("      - Cluster sizes:      from %d to %d (median %.0f)\n", min(cl_sizes_multi), max(cl_sizes_multi), median(cl_sizes_multi)))
  cat(sprintf("      - Singletons:         %d\n", sum(cl_sizes_multi == 1)))
}
cat(sprintf("\nGMYC analysis completed. Results in:\n   %s\n", OUT_DIR))
cat("   Delimitation_GMYC_single.csv\n")
cat("   Delimitation_GMYC_multi.csv\n")