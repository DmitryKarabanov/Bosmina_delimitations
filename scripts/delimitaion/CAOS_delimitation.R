# ==========================================================
# UNIVERSAL CAOS PIPELINE
# Comparison of 5 delimitation schemes: Ward (baseline) vs ASAP vs bGMYC vs KoT vs bPTP
# ==========================================================

DATA_DIR <- "C:/GEN/Bosmina/Delimitation/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

FASTA_FILE <- file.path(DATA_DIR, "_bosmina_all_GB_MEG_DJT.fasta")

# Delimitation schemes
WARD_CSV   <- file.path(DATA_DIR, "Delimitation_Ward_pair.csv")
ASAP_CSV   <- file.path(DATA_DIR, "Delimitation_ASAP_patri.csv")
BGMYC_CSV  <- file.path(DATA_DIR, "Delimitation_bGMYC_005.csv")
KOT_CSV    <- file.path(DATA_DIR, "Delimitation_KoT_t4.csv")
BPTP_CSV   <- file.path(DATA_DIR, "Delimitation_bPTP_p0.01.csv")

# Output files
OUT_DETAILS <- file.path(OUT_DIR, "CAOS_Universal_details.csv")
OUT_SUMMARY <- file.path(OUT_DIR, "CAOS_Universal_summary.csv")
OUT_CONFLICTS <- file.path(OUT_DIR, "CAOS_Universal_conflicts.csv")

# 1. Packages
pkgs <- c("ape", "dplyr", "tidyr")
lapply(pkgs, function(p) {
  if(!require(p, character.only = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  suppressMessages(library(p, character.only = TRUE))
})
cat("Packages loaded. Starting universal CAOS pipeline...\n")

# =====================================================================
# Data reading and synchronization
# =====================================================================
dna <- read.dna(FASTA_FILE, format = "fasta")

ward_df   <- read.table(WARD_CSV, sep = ";", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
asap_df   <- read.table(ASAP_CSV, sep = ";", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
bgmyc_df  <- read.table(BGMYC_CSV, sep = ";", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
kot_df    <- read.table(KOT_CSV, sep = ";", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
bptp_df   <- read.table(BPTP_CSV, sep = ";", header = TRUE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Merging all schemes
merged_df <- ward_df %>%
  inner_join(asap_df, by = "Sequence") %>%
  inner_join(bgmyc_df, by = "Sequence") %>%
  inner_join(kot_df, by = "Sequence") %>%
  inner_join(bptp_df, by = "Sequence")

colnames(merged_df) <- c("Sequence", "MOTU_Ward", "MOTU_ASAP", "MOTU_bGMYC", "MOTU_KoT", "MOTU_bPTP")

cat(sprintf("Loaded from FASTA: %d sequences\n", nrow(dna)))
cat(sprintf("After merging schemes: %d sequences\n", nrow(merged_df)))

# Outgroup filtering
OUTGROUP_PATTERN <- "Bosminopsis|Outgroup"
outgroup_idx <- grepl(OUTGROUP_PATTERN, rownames(dna), ignore.case = TRUE)
if(sum(outgroup_idx) > 0) {
  dna <- dna[!outgroup_idx, , drop = FALSE]
  merged_df <- merged_df[!grepl(OUTGROUP_PATTERN, merged_df$Sequence, ignore.case = TRUE), ]
  cat(sprintf("Removed outgroups: %d taxa\n", sum(outgroup_idx)))
}

# Strict synchronization
common_seqs <- intersect(rownames(dna), merged_df$Sequence)
dna <- dna[common_seqs, , drop = FALSE]
merged_df <- merged_df[match(common_seqs, merged_df$Sequence), ]

dna_matrix <- as.character(dna)
n_seqs <- nrow(dna_matrix)
n_sites <- ncol(dna_matrix)
cat(sprintf("Final set: %d taxa, %d sites\n", n_seqs, n_sites))

# =====================================================================
# FUNCTION: Search for simple diagnostic characters (Simple CA)
# =====================================================================
find_simple_CA <- function(motu_id, dna_matrix, clusters, n_sites) {
  in_motu_idx <- which(clusters == motu_id)
  out_motu_idx <- which(clusters != motu_id)
  
  if(length(in_motu_idx) == 0) return(list())
  
  diagnostic_sites <- list()
  for(site in 1:n_sites) {
    alleles_in <- dna_matrix[in_motu_idx, site]
    alleles_in <- alleles_in[!alleles_in %in% c("-", "n", "N", "?", "")]
    if(length(unique(alleles_in)) != 1) next
    
    diagnostic_base <- unique(alleles_in)
    alleles_out <- dna_matrix[out_motu_idx, site]
    alleles_out <- alleles_out[!alleles_out %in% c("-", "n", "N", "?", "")]
    
    if(diagnostic_base %in% alleles_out) next
    
    diagnostic_sites[[length(diagnostic_sites) + 1]] <- list(
      site = site, base = diagnostic_base,
      n_in = length(alleles_in), n_out = length(alleles_out)
    )
  }
  return(diagnostic_sites)
}

# =====================================================================
# SEARCH FOR DIAGNOSTIC SITES FOR THE BASELINE SCHEME (Ward)
# =====================================================================
cat("\nSearching for diagnostic sites for the baseline scheme (Ward)...\n")
unique_ward_clusters <- sort(unique(merged_df$MOTU_Ward))

results_list <- lapply(unique_ward_clusters, function(motu_id) {
  n_seqs_in <- sum(merged_df$MOTU_Ward == motu_id)
  diag_sites <- find_simple_CA(motu_id, dna_matrix, merged_df$MOTU_Ward, n_sites)
  
  # Progress for the first 10
  if(as.integer(motu_id) <= 10) {
    cat(sprintf("   Ward MOTU %s (%d seq): %d diagnostic sites\n",
                motu_id, n_seqs_in, length(diag_sites)))
  }
  
  list(
    motu_ward = motu_id,
    motu_asap = unique(merged_df$MOTU_ASAP[merged_df$MOTU_Ward == motu_id])[1],
    motu_bgmyc = unique(merged_df$MOTU_bGMYC[merged_df$MOTU_Ward == motu_id])[1],
    motu_kot = unique(merged_df$MOTU_KoT[merged_df$MOTU_Ward == motu_id])[1],
    motu_bptp = unique(merged_df$MOTU_bPTP[merged_df$MOTU_Ward == motu_id])[1],
    n_sequences = n_seqs_in,
    n_diagnostic_sites = length(diag_sites),
    sites = diag_sites
  )
})

# =====================================================================
# SUMMARY STATISTICS FOR ALL SCHEMES
# =====================================================================
cat("\nSummary statistics for all schemes:\n")
cat("------------------------------------------------------------\n")

summary_all <- list()
for(method in c("Ward", "ASAP", "bGMYC", "KoT", "bPTP")) {
  method_col <- paste0("motu_", tolower(method))
  if(method == "bGMYC") method_col <- "motu_bgmyc"
  if(method == "bPTP") method_col <- "motu_bptp"
  
  # Grouping by method MOTU
  method_summary <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      Method = method,
      MOTU = x[[method_col]],
      N_sequences = x$n_sequences,
      N_diag_sites = x$n_diagnostic_sites,
      stringsAsFactors = FALSE
    )
  }))
  
  # Aggregating: summing diagnostic sites for each MOTU
  method_agg <- method_summary %>%
    group_by(Method, MOTU) %>%
    summarise(
      Total_sequences = sum(N_sequences),
      Total_diag_sites = sum(N_diag_sites),
      Diagnosable = ifelse(sum(N_diag_sites) > 0, "Yes", "No"),
      .groups = 'drop'
    )
  
  n_total <- nrow(method_agg)
  n_diagnosable <- sum(method_agg$Diagnosable == "Yes")
  n_not_diagnosable <- n_total - n_diagnosable
  
  cat(sprintf("   %s: %d mOTUs | Diagnosable: %d (%.1f%%) | No CA: %d (%.1f%%)\n",
              method, n_total, n_diagnosable, 100 * n_diagnosable / n_total,
              n_not_diagnosable, 100 * n_not_diagnosable / n_total))
  
  summary_all[[method]] <- method_agg
}
cat("------------------------------------------------------------\n")

# Combining all summaries
summary_df <- do.call(rbind, summary_all)
write.table(summary_df, file = OUT_SUMMARY, row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Summary saved: %s\n", basename(OUT_SUMMARY)))

# =====================================================================
# CONFLICT ANALYSIS (Overlumping by conservative methods)
# =====================================================================
cat("\nCONFLICT ANALYSIS (when conservative methods 'lumped' valid lineages):\n")

base_df <- do.call(rbind, lapply(results_list, function(x) {
  data.frame(
    MOTU_Ward = x$motu_ward,
    MOTU_ASAP = x$motu_asap,
    MOTU_bGMYC = x$motu_bgmyc,
    MOTU_KoT = x$motu_kot,
    MOTU_bPTP = x$motu_bptp,
    N_sequences = x$n_sequences,
    N_diag_sites = x$n_diagnostic_sites,
    stringsAsFactors = FALSE
  )
}))

conflicts_list <- list()

# Checking each conservative method
for(method in c("ASAP", "bGMYC", "KoT")) {
  method_col <- paste0("MOTU_", method)
  if(method == "bGMYC") method_col <- "MOTU_bGMYC"
  
  conflict_analysis <- base_df %>%
    group_by(!!sym(method_col)) %>%
    summarise(
      n_ward_subgroups = n(),
      diagnosable_ward = sum(N_diag_sites > 0),
      ward_list = paste(MOTU_Ward, collapse = ", "),
      .groups = 'drop'
    ) %>%
    filter(n_ward_subgroups > 1 & diagnosable_ward > 1)
  
  if(nrow(conflict_analysis) > 0) {
    # FIX: Renaming the first column to a universal name
    names(conflict_analysis)[1] <- "Parent_MOTU"
    
    cat(sprintf("\n%s 'lumped' valid Ward lineages:\n", method))
    for(i in 1:nrow(conflict_analysis)) {
      cat(sprintf("   - %s MOTU %s contains %d Ward subgroups, of which %d have their own CA sites.\n", 
                  method,
                  conflict_analysis$Parent_MOTU[i], 
                  conflict_analysis$n_ward_subgroups[i],
                  conflict_analysis$diagnosable_ward[i]))
      cat(sprintf("     Ward composition: %s\n", conflict_analysis$ward_list[i]))
    }
    
    conflict_analysis$Method <- method
    conflicts_list[[method]] <- conflict_analysis
  } else {
    cat(sprintf("\n%s does not lump lineages that have unique CA sites.\n", method))
  }
}

# Saving conflicts
if(length(conflicts_list) > 0) {
  conflicts_df <- do.call(rbind, conflicts_list)
  # Reordering columns for clarity
  conflicts_df <- conflicts_df[, c("Method", "Parent_MOTU", "n_ward_subgroups", "diagnosable_ward", "ward_list")]
  write.table(conflicts_df, file = OUT_CONFLICTS, row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("\nConflicts saved: %s\n", basename(OUT_CONFLICTS)))
}

# =====================================================================
# DETAILED SITE TABLE
# =====================================================================
detail_list <- lapply(results_list, function(x) {
  if(x$n_diagnostic_sites > 0) {
    do.call(rbind, lapply(x$sites, function(s) {
      data.frame(
        MOTU_Ward = x$motu_ward,
        MOTU_ASAP = x$motu_asap,
        MOTU_bGMYC = x$motu_bgmyc,
        MOTU_KoT = x$motu_kot,
        MOTU_bPTP = x$motu_bptp,
        Site_Position = s$site,
        Diagnostic_Base = s$base,
        N_in = s$n_in,
        N_out = s$n_out,
        stringsAsFactors = FALSE
      )
    }))
  }
})

detail_df <- do.call(rbind, detail_list)
if(!is.null(detail_df)) {
  write.table(detail_df, file = OUT_DETAILS, row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("Detailed site table saved: %s\n", basename(OUT_DETAILS)))
} else {
  cat("Diagnostic sites not found. Detailed file not created.\n")
}

cat("\nUniversal CAOS pipeline successfully completed!\n")
cat("\nRESULTS:\n")
cat(sprintf("   %s - Summary for all methods\n", basename(OUT_SUMMARY)))
cat(sprintf("   %s - Detailed sites\n", basename(OUT_DETAILS)))
if(length(conflicts_list) > 0) {
  cat(sprintf("   %s - Conflicts (overlumping)\n", basename(OUT_CONFLICTS)))
}