# ==========================================================
# KoT — post-hoc validation of locMin scheme 
# ==========================================================

DATA_DIR <- "C:/GEN/Bosmina/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

FASTA_FILE     <- file.path(DATA_DIR, "_bosmina_all_GB_MEG_DJT.fasta")
START_CSV      <- file.path(OUT_DIR, "Delimitation_locMin_patri.csv")
THRESHOLDS     <- c(4, 6, 8)
OUT_CSV_LIST   <- file.path(OUT_DIR, paste0("Delimitation_KoT_t", THRESHOLDS, ".csv"))
OUT_MATRIX     <- file.path(OUT_DIR, "KoT_matrix_final.csv")
OUT_THETA      <- file.path(OUT_DIR, "KoT_theta_per_MOTU.csv")
OUT_HEATMAP    <- file.path(OUT_DIR, "KoT_heatmap.pdf")
OUT_DENDRO     <- file.path(OUT_DIR, "KoT_dendrogram.pdf")

# 1. Packages
pkgs <- c("ape")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting KoT (final version)...\n")

# =====================================================================
# VECTORIZED FUNCTIONS
# =====================================================================
watterson_theta <- function(dna_subset) {
  n <- nrow(dna_subset)
  if(n < 2) return(NA)
  a_n <- sum(1 / (1:(n - 1)))
  dna_matrix <- as.character(dna_subset)
  L <- ncol(dna_matrix)
  S <- sum(apply(dna_matrix, 2, function(col) {
    alleles <- col[!col %in% c("-", "n", "N", "?")]
    length(unique(alleles)) >= 2
  }))
  theta_W <- S / (a_n * L)
  return(theta_W)
}

count_S <- function(dna_subset) {
  dna_matrix <- as.character(dna_subset)
  S <- sum(apply(dna_matrix, 2, function(col) {
    alleles <- col[!col %in% c("-", "n", "N", "?")]
    length(unique(alleles)) >= 2
  }))
  return(S)
}

# =====================================================================
# Reading data
# =====================================================================
dna   <- read.dna(FASTA_FILE, format = "fasta")
start <- read.table(START_CSV, sep = ";", header = TRUE,
                    stringsAsFactors = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Loaded: %d sequences\n", nrow(dna)))
cat(sprintf("Starting scheme: %s (%d MOTUs)\n",
            basename(START_CSV), length(unique(start$MOTU_locMin))))

# Outgroup filtering (Bosminopsis)
OUTGROUP_PATTERN <- "Bosminopsis"
outgroup_idx <- grepl(OUTGROUP_PATTERN, rownames(dna), ignore.case = TRUE)
if(sum(outgroup_idx) > 0) {
  dna <- dna[!outgroup_idx, , drop = FALSE]
  cat(sprintf("Removed outgroups: %d taxa\n", sum(outgroup_idx)))
}

# Synchronization
dna_names   <- rownames(dna)
start_names <- start$Sequence
common <- intersect(dna_names, start_names)
dna   <- dna[common, , drop = FALSE]
start <- start[match(common, start$Sequence), ]
cat(sprintf("Taxa for analysis: %d\n", nrow(dna)))

clusters <- as.character(start$MOTU_locMin)
names(clusters) <- start$Sequence
unique_clusters <- sort(unique(clusters))
n_start <- length(unique_clusters)
cat(sprintf("Starting MOTUs (after filtering): %d\n", n_start))

L_alignment <- ncol(as.character(dna))
cat(sprintf("Alignment length: %d sites\n", L_alignment))

# =====================================================================
# CALCULATING theta_W FOR EACH STARTING MOTU
# =====================================================================
cat("\nCalculating theta_W (Watterson per site) for each MOTU...\n")
theta_per_cluster <- numeric(length(unique_clusters))
names(theta_per_cluster) <- unique_clusters
S_per_cluster <- numeric(length(unique_clusters))
names(S_per_cluster) <- unique_clusters
n_per_cluster <- numeric(length(unique_clusters))
names(n_per_cluster) <- unique_clusters

for(i in seq_along(unique_clusters)) {
  cl <- unique_clusters[i]
  idx <- which(clusters == cl)
  n_per_cluster[cl] <- length(idx)
  dna_subset <- dna[idx, , drop = FALSE]
  theta_per_cluster[cl] <- watterson_theta(dna_subset)
  S_per_cluster[cl] <- if(length(idx) >= 2) count_S(dna_subset) else 0
}

mean_theta <- mean(theta_per_cluster, na.rm = TRUE)
n_singletons_theta <- sum(is.na(theta_per_cluster))
theta_per_cluster[is.na(theta_per_cluster)] <- mean_theta
cat(sprintf("   Singletons (theta_W = NA -> mean theta): %d\n", n_singletons_theta))
cat(sprintf("   Mean theta_W: %.6f\n", mean_theta))
cat(sprintf("   Min theta_W: %.6f, Max theta_W: %.6f\n",
            min(theta_per_cluster), max(theta_per_cluster)))

theta_df <- data.frame(
  MOTU = unique_clusters,
  N_sequences = n_per_cluster,
  S_segregating_sites = S_per_cluster,
  Theta_Watterson = round(theta_per_cluster, 6),
  stringsAsFactors = FALSE
)
write.table(theta_df, file = OUT_THETA,
            row.names = FALSE, sep = ";", dec = ".",
            quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("theta_W table saved: %s\n", basename(OUT_THETA)))

# =====================================================================
# K MATRIX (inter-cluster K80 distances)
# =====================================================================
cat("\nCalculating K matrix (K80 inter-cluster distances)...\n")
dist_mat <- as.matrix(dist.dna(dna, model = "K80", pairwise.deletion = TRUE))
K_matrix <- matrix(0, nrow = n_start, ncol = n_start,
                   dimnames = list(unique_clusters, unique_clusters))

if(n_start >= 2) {
  for(i in 1:(n_start - 1)) {
    for(j in (i + 1):n_start) {
      idx_i <- which(clusters == unique_clusters[i])
      idx_j <- which(clusters == unique_clusters[j])
      K_val <- mean(dist_mat[idx_i, idx_j], na.rm = TRUE)
      K_matrix[i, j] <- K_val
      K_matrix[j, i] <- K_val
    }
  }
}

# =====================================================================
# K/theta MATRIX
# =====================================================================
cat("Calculating K/theta matrix...\n")
KoT_matrix <- matrix(NA, nrow = n_start, ncol = n_start,
                     dimnames = list(unique_clusters, unique_clusters))

if(n_start >= 2) {
  for(i in 1:(n_start - 1)) {
    for(j in (i + 1):n_start) {
      theta_avg <- (theta_per_cluster[i] + theta_per_cluster[j]) / 2
      if(theta_avg > 0 && !is.na(theta_avg) && K_matrix[i, j] > 0) {
        kot_val <- K_matrix[i, j] / theta_avg
        KoT_matrix[i, j] <- kot_val
        KoT_matrix[j, i] <- kot_val
      }
    }
  }
}

write.table(round(KoT_matrix, 2), file = OUT_MATRIX,
            sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("K/theta matrix saved: %s\n", basename(OUT_MATRIX)))

kot_values <- KoT_matrix[upper.tri(KoT_matrix)]
kot_values <- kot_values[!is.na(kot_values)]
cat(sprintf("\nK/theta statistics:\n"))
cat(sprintf("   Minimum:    %.2f\n", min(kot_values)))
cat(sprintf("   Median:     %.2f\n", median(kot_values)))
cat(sprintf("   Mean:       %.2f\n", mean(kot_values)))
cat(sprintf("   Maximum:    %.2f\n", max(kot_values)))

for(thr in THRESHOLDS) {
  n_below <- sum(kot_values < thr, na.rm = TRUE)
  pct <- 100 * n_below / length(kot_values)
  cat(sprintf("   Pairs with K/theta < %d: %d of %d (%.1f%%)\n",
              thr, n_below, length(kot_values), pct))
}

# =====================================================================
# COMPLETE-LINKAGE VIA hclust
# =====================================================================
cat("\nRunning complete-linkage clustering (hclust)...\n")
KoT_matrix_hc <- KoT_matrix
diag(KoT_matrix_hc) <- 0
max_kot <- max(KoT_matrix_hc, na.rm = TRUE)
KoT_matrix_hc[is.na(KoT_matrix_hc)] <- max_kot + 1
hc <- hclust(as.dist(KoT_matrix_hc), method = "complete")
cat("   Hierarchical tree built\n")

# =====================================================================
# DENDROGRAM WITH THRESHOLD LINES
# =====================================================================
cat("\nCreating dendrogram with threshold lines...\n")
tryCatch({
  pdf(file = OUT_DENDRO, width = 12, height = 8)
  par(mar = c(6, 4, 3, 8))
  plot(hc,
       main = "Hierarchical clustering of MOTU by K/theta (complete-linkage)",
       xlab = "MOTU (from locMin)", sub = "", cex = 0.8)
  
  thr_colors <- c("green3", "orange", "red")
  for(k in seq_along(THRESHOLDS)) {
    abline(h = THRESHOLDS[k], col = thr_colors[k], lty = 2, lwd = 2)
    text(x = n_start + 0.5, y = THRESHOLDS[k],
         labels = paste0("K/theta = ", THRESHOLDS[k],
                         " (p < ", c("0.05", "0.01", "0.001")[k], ")"),
         pos = 4, col = thr_colors[k], cex = 0.8, xpd = TRUE)
  }
  
  legend("topright", inset = c(-0.25, 0),
         legend = paste0("K/theta = ", THRESHOLDS, " -> ",
                         length(unique(cutree(hc, h = THRESHOLDS))), " MOTUs"),
         col = thr_colors, lty = 2, lwd = 2, cex = 0.8,
         bty = "n", xpd = TRUE)
  dev.off()
  cat(sprintf("   Dendrogram saved: %s\n", basename(OUT_DENDRO)))
}, error = function(e) {
  cat(sprintf("   Dendrogram error: %s\n", e$message))
})

# =====================================================================
# K/theta MATRIX HEATMAP
# =====================================================================
cat("\nCreating K/theta matrix heatmap...\n")
tryCatch({
  pdf(file = OUT_HEATMAP, width = 11, height = 10)
  n_colors <- 100
  pal <- colorRampPalette(c("darkblue", "white", "red"))(n_colors)
  
  KoT_viz <- KoT_matrix
  diag(KoT_viz) <- 0
  KoT_viz[is.na(KoT_viz)] <- max(KoT_viz, na.rm = TRUE) * 1.1
  
  if(requireNamespace("gplots", quietly = TRUE)) {
    gplots::heatmap.2(KoT_viz,
                      col = pal,
                      trace = "none",
                      dendrogram = "both",
                      margins = c(8, 8),
                      main = "K/theta matrix between MOTUs (locMin scheme)\nBlue = low (merge), Red = high (distinct)",
                      key.title = "K/theta",
                      cexRow = 0.8, cexCol = 0.8)
  } else {
    heatmap(KoT_viz, col = pal,
            main = "K/theta matrix between MOTUs\nBlue = low (merge), Red = high (distinct)",
            margins = c(8, 8))
  }
  dev.off()
  cat(sprintf("   Heatmap saved: %s\n", basename(OUT_HEATMAP)))
}, error = function(e) {
  cat(sprintf("   Heatmap error: %s\n", e$message))
  cat("   Hint: install gplots package for enhanced visualization\n")
})

# =====================================================================
# TREE CUTTING FOR THREE THRESHOLDS
# =====================================================================
cat("\nApplying K/theta thresholds...\n")
results_summary <- data.frame(
  Threshold = THRESHOLDS,
  Final_MOTU = rep(NA, length(THRESHOLDS)),
  Merges = rep(NA, length(THRESHOLDS)),
  Min_Size = rep(NA, length(THRESHOLDS)),
  Max_Size = rep(NA, length(THRESHOLDS)),
  Median_Size = rep(NA, length(THRESHOLDS)),
  Singletons = rep(NA, length(THRESHOLDS)),
  stringsAsFactors = FALSE
)

cluster_assignments <- list()

for(t_idx in seq_along(THRESHOLDS)) {
  thr <- THRESHOLDS[t_idx]
  cluster_idx <- cutree(hc, h = thr)
  names(cluster_idx) <- unique_clusters
  cluster_assignments[[t_idx]] <- cluster_idx
  
  new_clusters <- as.character(cluster_idx[clusters])
  names(new_clusters) <- names(clusters)
  n_final <- length(unique(new_clusters))
  
  cat(sprintf("\n   K/theta threshold = %d:\n", thr))
  cat(sprintf("      Starting MOTUs:     %d\n", n_start))
  cat(sprintf("      Final MOTUs:        %d\n", n_final))
  cat(sprintf("      Merges:             %d\n", n_start - n_final))
  
  cl_sizes <- table(new_clusters)
  cat(sprintf("      Cluster sizes:      from %d to %d (median %.0f)\n",
              min(cl_sizes), max(cl_sizes), median(cl_sizes)))
  cat(sprintf("      Singletons:         %d\n", sum(cl_sizes == 1)))
  
  results_summary$Final_MOTU[t_idx]  <- n_final
  results_summary$Merges[t_idx]      <- n_start - n_final
  results_summary$Min_Size[t_idx]    <- min(cl_sizes)
  results_summary$Max_Size[t_idx]    <- max(cl_sizes)
  results_summary$Median_Size[t_idx] <- median(cl_sizes)
  results_summary$Singletons[t_idx]  <- sum(cl_sizes == 1)
  
  unique_new <- unique(new_clusters)
  mapping <- setNames(seq_along(unique_new), unique_new)
  new_clusters_num <- as.integer(mapping[as.character(new_clusters)])
  names(new_clusters_num) <- names(new_clusters)
  
  kot_df <- data.frame(
    Sequence = names(new_clusters_num),
    MOTU_KoT = new_clusters_num,
    stringsAsFactors = FALSE
  )
  kot_df <- kot_df[order(kot_df$MOTU_KoT, kot_df$Sequence), ]
  rownames(kot_df) <- NULL
  
  write.table(kot_df, file = OUT_CSV_LIST[t_idx],
              row.names = FALSE, sep = ";", dec = ".",
              quote = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("      Saved: %s\n", basename(OUT_CSV_LIST[t_idx])))
  
  top5 <- sort(cl_sizes, decreasing = TRUE)[1:min(5, length(cl_sizes))]
  cat("      Top-5 largest:\n")
  for(motu in names(top5)) {
    cat(sprintf("         MOTU %s: %d seq\n", motu, top5[motu]))
  }
}

# =====================================================================
# MERGE MAP
# =====================================================================
cat("\nMerge map (which starting MOTUs were merged):\n")
cat("   (shows pairwise K/theta for merged MOTUs)\n")

for(t_idx in seq_along(THRESHOLDS)) {
  thr <- THRESHOLDS[t_idx]
  cluster_idx <- cluster_assignments[[t_idx]]
  
  cat(sprintf("\n   === K/theta threshold = %d (p < %s) ===\n",
              thr, c("0.05", "0.01", "0.001")[t_idx]))
  
  final_groups <- split(names(cluster_idx), cluster_idx)
  n_merged_groups <- sum(sapply(final_groups, length) > 1)
  n_singleton_groups <- sum(sapply(final_groups, length) == 1)
  
  cat(sprintf("      Clusters without merges: %d\n", n_singleton_groups))
  cat(sprintf("      Clusters after merges:   %d\n", n_merged_groups))
  cat("\n")
  
  for(grp_id in seq_along(final_groups)) {
    motus_in_group <- final_groups[[grp_id]]
    if(length(motus_in_group) == 1) {
      cat(sprintf("      KoT_%d <- [MOTU %s] (no merges)\n",
                  grp_id, motus_in_group))
    } else {
      cat(sprintf("      KoT_%d <- [%s]\n",
                  grp_id, paste(motus_in_group, collapse = " + ")))
      if(length(motus_in_group) >= 2) {
        for(i in 1:(length(motus_in_group) - 1)) {
          for(j in (i + 1):length(motus_in_group)) {
            mi <- motus_in_group[i]
            mj <- motus_in_group[j]
            idx_i <- which(unique_clusters == mi)
            idx_j <- which(unique_clusters == mj)
            kot <- KoT_matrix[idx_i, idx_j]
            cat(sprintf("            K/theta(%s, %s) = %.2f\n", mi, mj, kot))
          }
        }
      }
    }
  }
}

# =====================================================================
# COMPARATIVE SUMMARY
# =====================================================================
cat("\nComparative summary by thresholds:\n")
cat("   ================================================================\n")
cat(sprintf("   %-20s | %-8s | %-8s | %-8s\n", "Parameter",
            paste0("t=", THRESHOLDS[1]),
            paste0("t=", THRESHOLDS[2]),
            paste0("t=", THRESHOLDS[3])))
cat("   ---------------------------------------------------------------\n")
cat(sprintf("   %-20s | %8d | %8d | %8d\n", "Final MOTUs",
            results_summary$Final_MOTU[1], results_summary$Final_MOTU[2], results_summary$Final_MOTU[3]))
cat(sprintf("   %-20s | %8d | %8d | %8d\n", "Merges",
            results_summary$Merges[1], results_summary$Merges[2], results_summary$Merges[3]))
cat(sprintf("   %-20s | %8d | %8d | %8d\n", "Min size",
            results_summary$Min_Size[1], results_summary$Min_Size[2], results_summary$Min_Size[3]))
cat(sprintf("   %-20s | %8d | %8d | %8d\n", "Max size",
            results_summary$Max_Size[1], results_summary$Max_Size[2], results_summary$Max_Size[3]))
cat(sprintf("   %-20s | %8.0f | %8.0f | %8.0f\n", "Median size",
            results_summary$Median_Size[1], results_summary$Median_Size[2], results_summary$Median_Size[3]))
cat(sprintf("   %-20s | %8d | %8d | %8d\n", "Singletons",
            results_summary$Singletons[1], results_summary$Singletons[2], results_summary$Singletons[3]))
cat("   ================================================================\n")

cat("\nThreshold interpretation (according to Spori et al., 2021):\n")
cat("   K/theta = 4 (p < 0.05) — standard criterion\n")
cat("   K/theta = 6 (p < 0.01) — strict criterion\n")
cat("   K/theta = 8 (p < 0.001) — very strict criterion\n")
cat("\n   The HIGHER the threshold, the STRICTER the criterion 'this is one species',\n")
cat("   the MORE MOTUs are merged, the FEWER final MOTUs.\n")

cat(sprintf("\nKoT (final) completed. Results in:\n   %s\n", OUT_DIR))
cat("   Main tables:\n")
for(f in c(OUT_THETA, OUT_MATRIX, OUT_CSV_LIST)) {
  cat(sprintf("      %s\n", basename(f)))
}
cat("   Visualizations:\n")
cat(sprintf("      %s\n", basename(OUT_DENDRO)))
cat(sprintf("      %s\n", basename(OUT_HEATMAP)))