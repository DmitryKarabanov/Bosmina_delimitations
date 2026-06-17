# ==========================================================
# INTERACTIVE HEATMAP OF DELIMITATIONS + TREE
# Branch coloring by MOTU bGMYC (without legend)
# ==========================================================

# 1. Packages
pkgs  <- c("ape", "ggtree", "ggplot2", "dplyr", "tidyr", "purrr",
           "stringr", "plotly", "htmlwidgets", "RColorBrewer",
           "treeio", "tibble", "viridis")
for(p in pkgs) {
  if(!requireNamespace(p, quietly = TRUE)) {
    if(p %in% c("ggtree", "treeio")) {
      if(!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(p, ask = FALSE)
    } else {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
  suppressMessages(library(p, character.only = TRUE))
}
cat("All packages loaded\n")

# 2. Paths
DATA_DIR     <- "C:/GEN/Bosmina/"
DELIM_DIR    <- file.path(DATA_DIR, "Delimitation")
TREE_FILE    <- file.path(DATA_DIR, "_bosmina_all_GB_MEG_DJT.fasta.treefile")
BGMYC_FILE   <- file.path(DELIM_DIR, "Delimitation_bGMYC_005.csv")
OUT_HTML     <- file.path(DELIM_DIR, "Delimitation_heatmap_bgmyc_tree.html")
OUT_MATRIX   <- file.path(DELIM_DIR, "Delimitation_matrix_cleaned.csv")

# ==========================================================
# 3. READING AND CLEANING DELIMITATION FILES
# ==========================================================
cat("\nReading delimitation files...\n")
delim_files <- list.files(DELIM_DIR, pattern = "^Delimitation_.*\\.csv$", full.names = TRUE)
delim_files <- delim_files[!grepl("matrix_cleaned|congruence_summary|Agreement|bgmyc_tree", basename(delim_files))]
cat(sprintf("   Found files: %d\n", length(delim_files)))

read_and_clean  <- function(filepath) {
  df  <- read.table(filepath, sep = ";", header = TRUE,
                    stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  if (!"Sequence" %in% colnames(df)) return(NULL)
  
  method  <- basename(filepath) %>% str_remove("^Delimitation_") %>% str_remove("\\.csv$")
  motu_col  <- names(df)[2]
  
  outgroup_mask  <- grepl("Bosminopsis", df$Sequence, ignore.case = TRUE)
  if(sum(outgroup_mask) > 0) df  <- df[!outgroup_mask, ]
  
  df %>%
    select(Sequence, MOTU = all_of(motu_col)) %>%
    mutate(Method = method, MOTU = as.character(MOTU))
}

delim_list <- map(delim_files, read_and_clean) %>% compact()
delim_all  <- bind_rows(delim_list)
n_methods  <- length(unique(delim_all$Method))
cat(sprintf("   Methods: %d, Taxa: %d\n", n_methods, length(unique(delim_all$Sequence))))

# ==========================================================
# 4. TAXON x METHOD MATRIX
# ==========================================================
cat("\nBuilding matrix...\n")
delim_wide <- delim_all %>%
  pivot_wider(names_from = Method, values_from = MOTU) %>%
  as.data.frame()
rownames(delim_wide) <- delim_wide$Sequence
delim_wide$Sequence <- NULL

method_order  <- c(
  "ASAP_patri", "ASAP_pair", "ABGD_init_patri", "ABGD_init_pair",
  "ABGD_recur_patri", "ABGD_recur_pair", "locMin_patri", "locMin_pair",
  "Ward_patri", "Ward_pair", "Refined_patri", "Refined_pair",
  "GMYC_single", "GMYC_multi", "bGMYC_001", "bGMYC_005", "mPTP", "bPTP_p0.01",
  "KoT_t4", "KoT_t6", "KoT_t8"
)

available  <- intersect(method_order, colnames(delim_wide))
extra  <- setdiff(colnames(delim_wide), method_order)
final_method_order  <- c(available, extra)
delim_wide  <- delim_wide[, final_method_order]

write.table(delim_wide %>% rownames_to_column("Sequence"), file = OUT_MATRIX,
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")

# ==========================================================
# 5. IQ-TREE
# ==========================================================
cat("\nLoading tree...\n")
tree <- read.tree(TREE_FILE)
tree <- drop.tip(tree, tree$tip.label[grepl("Bosminopsis", tree$tip.label, ignore.case = TRUE)])

common_taxa <- intersect(tree$tip.label, rownames(delim_wide))
tree <- keep.tip(tree, common_taxa)
delim_wide <- delim_wide[common_taxa, , drop = FALSE]

p_tree <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tree_data <- p_tree$data

tips_data <- tree_data %>% filter(isTip) %>% arrange(y)
tip_order <- tips_data$label
delim_wide <- delim_wide[tip_order, , drop = FALSE]

n_tips <- length(tip_order)
n_methods <- ncol(delim_wide)

# ==========================================================
# 6. TREE COLORING BY bGMYC (O(N) ALGORITHM)
# ==========================================================
cat("\nLoading bGMYC and calculating branch colors...\n")
bgmyc_df <- read.table(BGMYC_FILE, sep = ";", header = TRUE,
                       stringsAsFactors = FALSE, fileEncoding = "UTF-8")
bgmyc_df <- bgmyc_df %>% filter(!grepl("Bosminopsis", Sequence, ignore.case = TRUE))

motu_col_bgmyc <- names(bgmyc_df)[2]
taxon_to_bgmyc <- setNames(as.character(bgmyc_df[[motu_col_bgmyc]]), bgmyc_df$Sequence)
unique_motus <- sort(unique(taxon_to_bgmyc))

# Palette for tree: viridis for MOTU + gray for Mixed
tree_motu_colors <- setNames(viridis::viridis(length(unique_motus), option = "plasma"), unique_motus)
tree_motu_colors["Mixed"] <- "#CCCCCC"

# O(N) bottom-up traversal
node_motu <- setNames(rep(NA, nrow(tree_data)), tree_data$node)
tip_mask <- tree_data$isTip

for(i in which(tip_mask)) {
  node_motu[as.character(tree_data$node[i])] <- taxon_to_bgmyc[tree_data$label[i]]
}

internal_nodes <- sort(tree_data$node[!tree_data$isTip], decreasing = TRUE)
for(node in internal_nodes) {
  children <- tree_data$node[tree_data$parent == node]
  children_motus <- node_motu[as.character(children)]
  children_motus <- children_motus[!is.na(children_motus)]
  
  if(length(children_motus) == 0) node_motu[as.character(node)] <- "Mixed"
  else if(length(unique(children_motus)) == 1) node_motu[as.character(node)] <- unique(children_motus)
  else node_motu[as.character(node)] <- "Mixed"
}

# ==========================================================
# 7. COLOR PALETTE FOR HEATMAP (SEPARATE!)
# ==========================================================
all_motus_heatmap <- sort(unique(unlist(lapply(delim_wide, unique))))
all_motus_heatmap <- all_motus_heatmap[!is.na(all_motus_heatmap)]
motu_palette_heatmap <- setNames(viridis::viridis(length(all_motus_heatmap), option = "plasma"), all_motus_heatmap)

# ==========================================================
# 8. METHOD LABELS (SHORTENED)
# ==========================================================
short_names  <- final_method_order %>%
  str_replace("_patri", " (pat)") %>% str_replace("_pair", " (K80)") %>%
  str_replace("_init", " init") %>% str_replace("_recur", " recur") %>%
  str_replace("_single", " single") %>% str_replace("_multi", " multi") %>%
  str_replace("_p0.01", " ")

# ==========================================================
# 9. PLOTLY: TREE (left part, colored by bGMYC)
# ==========================================================
cat("\nDrawing tree (colored by bGMYC, NO legend)...\n")
tree_edges <- tree_data %>% filter(!is.na(parent))
tree_fig <- plot_ly()
max_x_tree <- max(tree_data$x, na.rm = TRUE)

for(i in 1:nrow(tree_edges)) {
  row <- tree_edges[i, ]
  parent_row <- tree_data[tree_data$node == row$parent, ]
  if(nrow(parent_row) == 0) next
  
  child_motu <- node_motu[as.character(row$node)]
  branch_color <- ifelse(is.na(child_motu), "#CCCCCC", tree_motu_colors[child_motu])
  bl <- ifelse(is.na(row$branch.length), "N/A", round(row$branch.length, 5))
  
  if(row$isTip && !is.na(row$label)) {
    hover_txt  <- paste0("<b>Taxon: </b>", row$label,
                         "<br><b>bGMYC MOTU: </b>", child_motu,
                         "<br><b>Branch length: </b>", bl)
  } else {
    desc_count  <- sum(tree_data$parent == row$node, na.rm = TRUE)
    hover_txt  <- paste0("<b>Node: </b>", row$node,
                         "<br><b>bGMYC MOTU: </b>", child_motu,
                         "<br><b>Children: </b>", desc_count,
                         "<br><b>Branch length: </b>", bl)
  }
  
  # Thick branches: width = 5.4
  tree_fig <- tree_fig %>% add_segments(
    x = parent_row$x, xend = row$x, y = row$y, yend = row$y,
    line = list(color = branch_color, width = 5.4),
    hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE
  )
  
  tree_fig <- tree_fig %>% add_segments(
    x = parent_row$x, xend = parent_row$x, y = parent_row$y, yend = row$y,
    line = list(color = branch_color, width = 5.4),
    hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE
  )
}

# Dots at branch ends (colored by bGMYC MOTU)
tree_fig  <- tree_fig %>% add_markers(
  data = tips_data, x = ~x, y = ~y,
  marker = list(size = 3,
                color = sapply(tips_data$label, function(l) {
                  m  <- taxon_to_bgmyc[l]
                  ifelse(is.na(m), "#CCCCCC", tree_motu_colors[m])
                }),
                line = list(width = 0)),
  hovertext = ~paste0("<b>", label, "</b><br><b>bGMYC: </b>",
                      sapply(label, function(l) taxon_to_bgmyc[l])),
  hoverinfo = "text", showlegend = FALSE
)

# Large sparse labels
label_step <- max(1, ceiling(n_tips / 60))
tips_labeled <- tips_data %>%
  mutate(row_num = row_number()) %>%
  filter(row_num %% label_step == 1)

cat(sprintf("   Showing %d labels out of %d (every %d-th)\n",
            nrow(tips_labeled), n_tips, label_step))

label_offset <- max_x_tree * 0.25
tree_fig <- tree_fig %>% add_text(
  data = tips_labeled,
  x = ~x + label_offset,
  y = ~y,
  text = ~label,
  textposition = "middle left",
  textfont = list(size = 18, color = "black", family = "monospace"),
  showlegend = FALSE,
  hoverinfo = "skip"
)

tree_fig <- tree_fig %>% layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(min(tips_data$y) - 1, max(tips_data$y) + 1))
)

# ==========================================================
# 10. PLOTLY: HEATMAP (right part)
# ==========================================================
cat("Drawing heatmap...\n")
label_space <- max_x_tree * 0.4
heatmap_offset <- max_x_tree + label_space
cell_width <- max_x_tree * 0.025
cell_height <- 0.8

shapes_list <- list()
hover_df <- data.frame(x = numeric(), y = numeric(), text = character(), stringsAsFactors = FALSE)

for(j in 1:n_methods) {
  method_name <- colnames(delim_wide)[j]
  x_left <- heatmap_offset + (j - 1) * cell_width
  x_right <- heatmap_offset + j * cell_width
  x_mid <- (x_left + x_right) / 2
  
  method_col <- as.character(delim_wide[[j]])
  
  for(i in 1:n_tips) {
    taxon <- rownames(delim_wide)[i]
    motu_val <- method_col[i]
    tip_row <- tips_data[tips_data$label == taxon, ]
    y_pos <- tip_row$y
    
    cell_color <- if(is.na(motu_val)) "#F0F0F0" else motu_palette_heatmap[motu_val]
    
    shapes_list[[length(shapes_list) + 1]]  <- list(
      type = "rect", x0 = x_left, x1 = x_right,
      y0 = y_pos - cell_height/2, y1 = y_pos + cell_height/2,
      fillcolor = cell_color, line = list(color = "white", width = 0.5)
    )
    
    hover_txt  <- sprintf("<b>Taxon: </b>%s<br><b>Method: </b>%s<br><b>MOTU: </b>%s",
                          taxon, method_name, ifelse(is.na(motu_val), "N/A", motu_val))
    hover_df  <- rbind(hover_df, data.frame(x = x_mid, y = y_pos, text = hover_txt))
  }
  
  runs  <- rle(method_col)
  y_start_idx  <- 1
  for(k in seq_along(runs$lengths)) {
    run_len  <- runs$lengths[k]
    run_val  <- runs$values[k]
    y_end_idx  <- y_start_idx + run_len - 1
    
    if(!is.na(run_val)) {
      taxa_in_run  <- rownames(delim_wide)[y_start_idx:y_end_idx]
      y_positions  <- tips_data$y[match(taxa_in_run, tips_data$label)]
      y_min  <- min(y_positions) - cell_height/2
      y_max  <- max(y_positions) + cell_height/2
      
      shapes_list[[length(shapes_list) + 1]]  <- list(
        type = "rect", x0 = x_left, x1 = x_right, y0 = y_min, y1 = y_max,
        fillcolor = "rgba(0,0,0,0)", line = list(color = "black", width = 2), layer = "above"
      )
    }
    y_start_idx  <- y_end_idx + 1
  }
}

matrix_fig <- plotly_empty() %>% layout(shapes = shapes_list)
matrix_fig <- matrix_fig %>% add_trace(
  type = "scatter", mode = "markers", data = hover_df, x = ~x, y = ~y,
  hovertext = ~text, hoverinfo = "text",
  marker = list(size = 18, color = "rgba(0,0,0,0)", line = list(width = 0)),
  showlegend = FALSE
)

for(j in 1:n_methods) {
  matrix_fig <- matrix_fig %>% add_annotations(
    x = heatmap_offset + (j - 0.5) * cell_width, y = max(tips_data$y) + 3,
    text = short_names[j], textangle = 90,
    font = list(size = 9, color = "black"), showarrow = FALSE
  )
}

matrix_fig <- matrix_fig %>% layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(min(tips_data$y) - 1, max(tips_data$y) + 1))
)

# ==========================================================
# 11. COMBINING VIA SUBPLOT
# ==========================================================
cat("Combining tree and heatmap...\n")
combined <- subplot(tree_fig, matrix_fig, nrows = 1, widths = c(0.45, 0.55),
                    shareY = TRUE, titleX = FALSE, titleY = FALSE)

x_max <- heatmap_offset + n_methods * cell_width + 0.1
combined  <- combined %>% layout(
  title = list(
    text = paste0("<b>Integrative Delimitation of <i>Bosmina</i>: Phylogeny + Methods</b><br>",
                  "<span style='font-size:12px'>", n_tips, " taxa x ", n_methods, " methods | ",
                  "Branch colors = bGMYC MOTUs (gray = mixed) | Black boxes = MOTU boundaries</span>"),
    x = 0.5
  ),
  xaxis = list(range = c(-0.02, x_max), showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
  plot_bgcolor = "white", hovermode = "closest",
  margin = list(l = 10, r = 20, t = 100, b = 10)
)

# ==========================================================
# 12. SAVING
# ==========================================================
htmlwidgets::saveWidget(combined, OUT_HTML, selfcontained = TRUE, title = "Tree and Delimitation")

cat("\n=================================================\n")
cat("  DONE!\n")
cat("=================================================\n")
cat(sprintf("  File: %s\n", OUT_HTML))
cat("\n  Features of this version:\n")
cat("     - Branches colored by consensus bGMYC MOTUs\n")
cat("     - Gray branches = nodes uniting different MOTUs (paraphyly)\n")
cat("     - Legend removed for plot clarity\n")
cat("     - Hover on branch shows bGMYC MOTU ID\n")
cat("     - Use Autoscale for optimal size\n")
cat("=================================================\n")