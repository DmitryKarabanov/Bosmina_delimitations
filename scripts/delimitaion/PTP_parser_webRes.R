# ==========================================================
# bPTP RESULTS PARSER 
# ==========================================================

# 1. Directory setup
DATA_DIR <- "C:/GEN/Bosmina/Delimitation/"
OUT_DIR  <- "C:/GEN/Bosmina/Delimitation/"
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# File search (several possible names)
BPTP_CANDIDATES <- c("bPTP_results.txt", "bptp_results_p0.01.txt",
                     "bptp_results.txt", "bPTP_p0.01.txt",
                     "bPTP_best_partition.txt")
BPTP_FILE <- NULL
for(f in BPTP_CANDIDATES) {
  path <- file.path(DATA_DIR, f)
  if(file.exists(path)) {
    BPTP_FILE <- path
    break
  }
}

if(is.null(BPTP_FILE)) {
  cat("bPTP file not found in standard locations. Please select it manually:\n")
  BPTP_FILE <- file.choose()
}

OUT_CSV <- file.path(OUT_DIR, "Delimitation_bPTP_p0.01.csv")
cat(sprintf("Reading file: %s\n", BPTP_FILE))

# 2. Reading file
lines_raw <- readLines(BPTP_FILE, encoding = "UTF-8", warn = FALSE)
lines <- trimws(gsub("^\uFEFF", "", lines_raw))
cat(sprintf("Total lines in file: %d\n", length(lines)))

# 3. Grouping: splitting by empty lines
# Each group of consecutive non-empty lines = one MOTU
groups <- list()
current_group <- c()

for(line in lines) {
  if(nchar(line) == 0) {
    # Empty line — save current group (if non-empty)
    if(length(current_group) > 0) {
      groups <- c(groups, list(current_group))
      current_group <- c()
    }
  } else {
    # Non-empty line — add to current group
    current_group <- c(current_group, line)
  }
}

# Don't forget the last group
if(length(current_group) > 0) {
  groups <- c(groups, list(current_group))
}
cat(sprintf("Found groups (MOTUs): %d\n", length(groups)))

# 4. Forming the table
records <- list()
for(motu_id in seq_along(groups)) {
  seq_names <- groups[[motu_id]]
  for(seq_name in seq_names) {
    records[[seq_name]] <- list(
      MOTU   = motu_id,
      Size   = length(seq_names)
    )
  }
}
cat(sprintf("Parsed sequences: %d\n", length(records)))

# 5. Creating data.frame (ONLY 2 columns)
sequences <- names(records)
motus     <- sapply(records, function(x) as.integer(x$MOTU))

bptp_df <- data.frame(
  Sequence  = as.character(sequences),
  MOTU_bPTP = as.integer(motus),
  stringsAsFactors = FALSE
)

# 6. Sorting: by MOTU, then by Sequence
bptp_df <- bptp_df[order(bptp_df$MOTU_bPTP, bptp_df$Sequence), ]
rownames(bptp_df) <- NULL

# 7. Saving (separator ";", UTF-8, no quotes)
write.table(bptp_df, file = OUT_CSV,
            row.names = FALSE, sep = ";", dec = ".",
            quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Table saved: %s\n", basename(OUT_CSV)))

# 8. Detailed summary
n_total <- nrow(bptp_df)
n_motu  <- length(groups)
cl_sizes <- table(bptp_df$MOTU_bPTP)

cat("\nSummary for bPTP (p=0.01):\n")
cat(sprintf("   - Taxa analyzed:          %d\n", n_total))
cat(sprintf("   - Total MOTUs:            %d\n", n_motu))
cat(sprintf("   - Cluster sizes:          from %d to %d (median %.0f)\n",
            min(cl_sizes), max(cl_sizes), median(cl_sizes)))
cat(sprintf("   - Singletons (1 seq):     %d (%.1f%%)\n",
            sum(cl_sizes == 1), 100 * sum(cl_sizes == 1) / n_motu))

# 9. Largest clusters (top-10)
cat("\nTop-10 largest MOTUs:\n")
top10 <- sort(cl_sizes, decreasing = TRUE)[1:min(10, length(cl_sizes))]
for(motu in names(top10)) {
  cat(sprintf("   - MOTU %2s: %3d seq\n", motu, top10[motu]))
}

# 10. Singletons
singletons <- as.integer(names(cl_sizes[cl_sizes == 1]))
if(length(singletons) > 0) {
  cat(sprintf("\nSingletons (%d):\n", length(singletons)))
  for(s in singletons) {
    seq_name <- bptp_df$Sequence[bptp_df$MOTU_bPTP == s][1]
    cat(sprintf("   - MOTU %d: %s\n", s, seq_name))
  }
}

# 11. Size distribution
cat("\nMOTU size distribution:\n")
size_dist <- table(cl_sizes)
for(sz in as.integer(names(size_dist))) {
  cat(sprintf("   - Size %2d: %2d clusters\n", sz, size_dist[as.character(sz)]))
}

cat(sprintf("\nParsing completed. Result: %s\n", OUT_CSV))