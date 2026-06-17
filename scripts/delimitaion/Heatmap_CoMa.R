# ==========================================================
# INTERACTIVE AGREEMENT MATRIX + COLORED TREE (v6)
# NO legend, font=18, dots=2, lines=5.4, labels on the right
# ==========================================================

# 1. Packages
pkgs  <- c("ape", "dplyr", "tidyr", "purrr", "stringr",
           "ggtree", "treeio", "plotly", "htmlwidgets")
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
cat("Packages loaded\n")

# 2. Paths
DATA_DIR  <- "C:/GEN/Bosmina/"
DELIM_DIR <- file.path(DATA_DIR, "Delimitation")
TREE_FILE <- file.path(DATA_DIR, "_bosmina_all_GB_MEG_DJT.fasta.treefile")
OUT_HTML  <- file.path(DELIM_DIR, "Agreement_Matrix.html")

# ==========================================================
# 3. READING DELIMITATION FILES
# ==========================================================
cat("\nReading files...\n")
delim_files <- list.files(DELIM_DIR, pattern = "^Delimitation_.*\\.csv$", full.names = TRUE)
delim_files <- delim_files[!grepl("matrix_cleaned|congruence_summary|Agreement", basename(delim_files))]
cat(sprintf("   Found files: %d\n", length(delim_files)))

read_and_clean <- function(filepath) {
  df  <- read.table(filepath, sep = ";", header = TRUE,
                    stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  if (!"Sequence" %in% colnames(df)) return(NULL)
  
  method  <- basename(filepath) %>% str_remove("^Delimitation_") %>% str_remove("\\.csv$")
  df  <- df %>% filter(!grepl("Bosminopsis", Sequence, ignore.case = TRUE))
  motu_col  <- names(df)[2]
  
  df %>%
    select(Sequence, MOTU = all_of(motu_col)) %>%
    mutate(Method = method, MOTU = as.character(MOTU))
}

delim_list <- map(delim_files, read_and_clean) %>% compact()
all_data   <- bind_rows(delim_list)
n_methods  <- length(unique(all_data$Method))
cat(sprintf("   Methods: %d, Taxa: %d\n", n_methods, length(unique(all_data$Sequence))))

# ==========================================================
# 4. N x N AGREEMENT MATRIX
# ==========================================================
cat("\nCalculating agreement matrix...\n")
wide_data <- all_data %>%
  pivot_wider(names_from = Method, values_from = MOTU) %>%
  as.data.frame()
rownames(wide_data) <- wide_data$Sequence
wide_data$Sequence <- NULL

taxa   <- rownames(wide_data)
n_taxa <- length(taxa)

agreement_matrix <- matrix(0, nrow = n_taxa, ncol = n_taxa,
                           dimnames = list(taxa, taxa))

for(i in 1:ncol(wide_data)) {
  motu_vec <- as.character(wide_data[[i]])
  bin_mat  <- outer(motu_vec, motu_vec,
                    function(a, b) ifelse(is.na(a) | is.na(b), 0, a == b))
  agreement_matrix <- agreement_matrix + bin_mat
}
agreement_matrix <- agreement_matrix / n_methods
diag(agreement_matrix) <- 1
cat(sprintf("   Matrix %dx%d ready\n", n_taxa, n_taxa))

# ==========================================================
# 5. TREE + OPTIMIZED CLADE ASSIGNMENT
# ==========================================================
cat("\nLoading tree...\n")
tree <- read.tree(TREE_FILE)
tree <- drop.tip(tree, tree$tip.label[grepl("Bosminopsis", tree$tip.label, ignore.case = TRUE)])

common_taxa <- intersect(taxa, tree$tip.label)
tree <- keep.tip(tree, common_taxa)

p_tree <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tree_data <- p_tree$data

tips_data <- tree_data %>% filter(isTip) %>% arrange(y)
tip_order <- tips_data$label
agreement_matrix <- agreement_matrix[tip_order, tip_order]
n_tips <- length(tip_order)

# Clade assignment function
assign_clade <- function(name) {
  if(grepl("Eubosmina|tubicen|coregoni|longispina|longirostris|liederi|fatalis|freyi|meridionalis", name, ignore.case = TRUE)) return("Eubosmina")
  if(grepl("Liederobosmina|huaronensis|chilense|tanakai", name, ignore.case = TRUE)) return("Liederobosmina")
  if(grepl("Lunobosmina|oriens", name, ignore.case = TRUE)) return("Lunobosmina")
  if(grepl("korineki|Colombia", name, ignore.case = TRUE)) return("Colombian_clade")
  return("Other")
}

clade_colors <- c(
  "Eubosmina"       = "#E41A1C",
  "Liederobosmina"  = "#377EB8",
  "Lunobosmina"     = "#4DAF4A",
  "Colombian_clade" = "#984EA3",
  "Other"           = "#FF7F00",
  "Mixed"           = "#999999"
)

tip_clades <- setNames(sapply(tips_data$label, assign_clade), tips_data$label)

# O(N) clade assignment
cat("Assigning clades (O(N) algorithm)...\n")
node_clade <- setNames(rep(NA, nrow(tree_data)), tree_data$node)
tip_mask   <- tree_data$isTip
tip_labels <- tree_data$label[tip_mask]
tip_nodes  <- tree_data$node[tip_mask]

for(i in seq_along(tip_nodes)) {
  node_clade[as.character(tip_nodes[i])] <- assign_clade(tip_labels[i])
}

internal_nodes <- tree_data$node[!tree_data$isTip]
internal_nodes <- sort(internal_nodes, decreasing = TRUE)

for(node in internal_nodes) {
  children <- tree_data$node[tree_data$parent == node]
  if(length(children) == 0) {
    node_clade[as.character(node)] <- "Mixed"
    next
  }
  children_clades <- node_clade[as.character(children)]
  children_clades <- children_clades[!is.na(children_clades)]
  
  if(length(children_clades) == 0) {
    node_clade[as.character(node)] <- "Mixed"
  } else if(length(unique(children_clades)) == 1) {
    node_clade[as.character(node)] <- unique(children_clades)
  } else {
    node_clade[as.character(node)] <- "Mixed"
  }
}

# ==========================================================
# 6. PLOTLY: TREE (left half)
# ==========================================================
cat("\nDrawing colored tree...\n")
tree_edges <- tree_data %>% filter(!is.na(parent))
tree_fig   <- plot_ly()
max_x_tree <- max(tree_data$x, na.rm = TRUE)

# Branches with color
for(i in 1:nrow(tree_edges)) {
  row <- tree_edges[i, ]
  parent_row <- tree_data[tree_data$node == row$parent, ]
  if(nrow(parent_row) == 0) next
  
  child_clade  <- node_clade[as.character(row$node)]
  branch_color <- ifelse(is.na(child_clade), "#999999", clade_colors[child_clade])
  bl <- ifelse(is.na(row$branch.length), "N/A", round(row$branch.length, 5))
  
  if(row$isTip && !is.na(row$label)) {
    hover_txt <- paste0(
      "<b>Taxon: </b>", row$label, "<br>",
      "<b>Clade: </b>", child_clade, "<br>",
      "<b>Branch: </b>", bl)
  } else {
    desc_count <- sum(tree_data$parent == row$node, na.rm = TRUE)
    hover_txt <- paste0(
      "<b>Node: </b>", row$node, "<br>",
      "<b>Clade: </b>", child_clade, "<br>",
      "<b>Direct children: </b>", desc_count, "<br>",
      "<b>Branch: </b>", bl)
  }
  
  tree_fig <- tree_fig %>% add_segments(
    x = parent_row$x, xend = row$x, y = row$y, yend = row$y,
    line = list(color = branch_color, width = 6),
    hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE
  )
  tree_fig <- tree_fig %>% add_segments(
    x = parent_row$x, xend = parent_row$x, y = parent_row$y, yend = row$y,
    line = list(color = branch_color, width = 6),
    hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE
  )
}

# Dots at branch ends
tree_fig <- tree_fig %>% add_markers(
  data = tips_data, x = ~x, y = ~y,
  marker = list(size = 2, color = "black", line = list(width = 0)),
  hovertext = ~paste0("<b>", label, "</b><br><b>Clade: </b>",
                      sapply(label, function(x) tip_clades[x])),
  hoverinfo = "text", showlegend = FALSE
)

# Labels on the right
label_step <- max(1, ceiling(n_tips / 60))
tips_labeled <- tips_data %>%
  mutate(row_num = row_number()) %>%
  filter(row_num %% label_step == 1)

label_offset <- max_x_tree * 0.25
cat(sprintf("   Showing %d labels out of %d (every %d-th)\n",
            nrow(tips_labeled), n_tips, label_step))

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
# 7. PLOTLY: MATRIX (right half, NO bottom labels)
# ==========================================================
cat("Drawing agreement matrix...\n")
hover_text <- matrix(
  paste0("<b>Taxon 1: </b>", rownames(agreement_matrix)[row(agreement_matrix)],
         "<br><b>Taxon 2: </b>", colnames(agreement_matrix)[col(agreement_matrix)],
         "<br><b>Agreement: </b>", sprintf("%.2f", agreement_matrix),
         " (", round(agreement_matrix * n_methods, 0), "/", n_methods, " methods)"),
  nrow = n_tips
)

matrix_fig <- plot_ly(
  z = agreement_matrix, x = 1:n_tips, y = 1:n_tips, type = "heatmap",
  colorscale = list(
    c(0,    "#F0F0F0"),
    c(0.3,  "#A6CEE3"),
    c(0.7,  "#1F78B4"),
    c(1,    "#08306B")
  ),
  zmin = 0, zmax = 1,
  text = hover_text, hoverinfo = "text", showscale = TRUE,
  colorbar = list(
    title = "Agreement",
    tickvals = c(0, 0.5, 1),
    ticktext = c("0%", "50%", "100%"),
    len = 0.6
  )
)

matrix_fig <- matrix_fig %>% layout(
  xaxis = list(title = "", showticklabels = FALSE, showgrid = FALSE, zeroline = FALSE, range = c(0.5, n_tips + 0.5)),
  yaxis = list(title = "", showticklabels = FALSE, autorange = "reversed", showgrid = FALSE, zeroline = FALSE, range = c(0.5, n_tips + 0.5))
)

# ==========================================================
# 8. COMBINING VIA SUBPLOT (NO LEGEND!)
# ==========================================================
cat("Combining tree and matrix...\n")
combined <- subplot(
  tree_fig, matrix_fig, nrows = 1, widths = c(0.50, 0.50),
  shareY = TRUE, titleX = FALSE, titleY = FALSE
)

combined <- combined %>% layout(
  title = list(
    text = paste0(
      "<b>Integrative Delimitation: Agreement Matrix + Clade-colored Phylogeny</b><br>",
      "<span style='font-size:12px'>", n_tips, " taxa x ", n_methods, " methods | ",
      "Hover branches for details | Colors = major clades</span>"
    ),
    x = 0.5
  ),
  margin = list(l = 20, r = 80, t = 100, b = 30)
)

# ==========================================================
# 9. SAVING
# ==========================================================
htmlwidgets::saveWidget(combined, OUT_HTML, selfcontained = TRUE, title = "Agreement Matrix + Colored Tree v6")

cat("\n=================================================\n")
cat("  DONE!\n")
cat("=================================================\n")