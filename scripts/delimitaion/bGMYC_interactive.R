#═══════════════════════════════════════════════════════════════
# bGMYC4 v4.1.0: Full interactive pipeline
#═══════════════════════════════════════════════════════════════

library(ape)
library(future)
library(future.apply)
library(mcmcse)
library(plotly)
library(base64enc)
library(htmlwidgets)
library(treeio)
library(dplyr)
library(ggtree)
library(bGMYC4)

cat(" All critical dependencies loaded\n")
cat("===  bGMYC4 Interactive Analysis ===\n\n")

# 1. LOAD TREES  
consensus_path <- readline(prompt = " Path to consensus tree (tree): ")
posterior_path <- readline(prompt = " Path to posterior tree set (trees): ")

if (!file.exists(consensus_path)) stop(" Consensus tree not found!")
if (!file.exists(posterior_path)) stop(" Trees file not found!")

DELIM_DIR <- dirname(consensus_path)
html_path <- file.path(DELIM_DIR, "bGMYC_interactive_heatmap.html")
cat(sprintf(" Results will be saved to: %s\n", DELIM_DIR))

# Use treeio::read.beast() to extract posterior from Nexus annotations
tree_beast <- read.beast(consensus_path)
tree <- tree_beast@phylo  # Extract phylo object for bGMYC

# For posterior_trees use standard read.nexus
all_trees <- read.nexus(posterior_path)
class(all_trees) <- "multiPhylo"

cat(sprintf(" Loaded: 1 consensus + %d trees in posterior\n", length(all_trees)))

# 2. SELECT ANALYSIS MODE 
cat("\n Select analysis mode:\n")
cat("   1. Consensus tree only (fast, ignores topological uncertainty)\n")
cat("   2. Consensus + multiple trees (accounts for phylogenetic uncertainty)\n")

mode_str <- readline("   Enter 1 or 2 (default 2): ")
analysis_mode <- ifelse(nchar(trimws(mode_str)) == 0 || is.na(as.integer(mode_str)), 2, as.integer(mode_str))

if (!analysis_mode %in% c(1, 2)) {
  cat(" Invalid choice. Setting mode 2.\n")
  analysis_mode <- 2
}

# 3. SUBSAMPLE TREES (ONLY FOR MODE 2)
if (analysis_mode == 2) {
  n_total <- length(all_trees)
  burnin_idx <- floor(n_total * 0.10)
  
  n_str <- readline(prompt = sprintf("\n How many trees to sample? (default 10, available: %d): ", n_total - burnin_idx))
  n_sample <- ifelse(nchar(n_str) == 0 || is.na(as.numeric(n_str)), 10, as.integer(n_str))
  
  if (n_sample > (n_total - burnin_idx)) n_sample <- n_total - burnin_idx
  
  set.seed(42)
  trees_sample <- all_trees[sample((burnin_idx + 1):n_total, n_sample)]
  class(trees_sample) <- "multiPhylo"
  
  cat(sprintf(" Sampled %d random trees (indices %d–%d)\n", n_sample, burnin_idx + 1, n_total))
} else {
  cat("\n Mode 1 selected: analysis on consensus tree only.\n")
  trees_sample <- NULL
}

# 4. CHECK ULTRAMETRICITY
fix_ultrametric <- function(tr) {
  if (!is.ultrametric(tr)) {
    tr$edge.length <- round(tr$edge.length, 8)
    if (!is.ultrametric(tr)) stop(" Tree remains non-ultrametric after rounding branch lengths.\n")
  }
  tr
}

cat(" Checking ultrametricity...\n")
tree <- fix_ultrametric(tree)

if (analysis_mode == 2) {
  trees_sample <- lapply(trees_sample, fix_ultrametric)
  class(trees_sample) <- "multiPhylo"
}

cat(" Trees are ultrametric.\n")

# 5. OUTGROUP (WITH PATTERNS)
og_prompt <- readline(prompt = "\n Remove outgroups? (y/n, default n): ")

if (tolower(og_prompt) == "y") {
  og_str <- readline(prompt = "   Enter outgroup names (comma-separated): ")
  og_patterns <- trimws(unlist(strsplit(og_str, "[,;]+")))
  og_patterns <- og_patterns[nchar(og_patterns) > 0]
  
  if (length(og_patterns) > 0) {
    matching_tips <- c()
    unmatched_patterns <- c()
    
    for(p in og_patterns) {
      pattern <- paste0("(^|[^A-Za-z0-9])", p, "($|[^A-Za-z0-9])")
      found <- tree$tip.label[grepl(pattern, tree$tip.label, ignore.case = TRUE, perl = TRUE)]
      
      if (length(found) > 0) {
        matching_tips <- c(matching_tips, found)
      } else {
        unmatched_patterns <- c(unmatched_patterns, p)
      }
    }
    
    matching_tips <- unique(matching_tips)
    
    if (length(unmatched_patterns) > 0) {
      cat(sprintf("  Not found in tree: %s\n", paste(unmatched_patterns, collapse = ", ")))
    }
    
    if (length(matching_tips) > 0) {
      tree <- drop.tip(tree, matching_tips)
      
      if (exists("analysis_mode") && exists("trees_sample")) {
        if (analysis_mode == 2) {
          trees_sample <- lapply(trees_sample, drop.tip, tip = matching_tips)
          class(trees_sample) <- "multiPhylo"
        }
      }
      
      cat(sprintf(" Removed %d taxa:\n", length(matching_tips)))
      if (length(matching_tips) > 10) {
        cat(sprintf("   %s ... and %d more\n", 
                    paste(head(matching_tips, 10), collapse = ", "), 
                    length(matching_tips) - 10))
      } else {
        cat(sprintf("   %s\n", paste(matching_tips, collapse = ", ")))
      }
    } else {
      cat("  No taxa matched the provided names. Skipping.\n")
    }
  }
}

ntips <- length(tree$tip.label)
cat(sprintf(" Final size: %d taxa\n\n", ntips))

# 6. PARAMETER INPUT BLOCK 
input_scalar <- function(label, default_val) {
  val <- readline(sprintf("  %s (default: %s): ", label, default_val))
  if (nchar(trimws(val)) == 0) return(default_val)
  num <- suppressWarnings(as.numeric(val))
  if (is.na(num)) {
    cat("  Non‑numeric input. Using default.\n")
    return(default_val)
  }
  return(num)
}

input_vector <- function(label, default_vec) {
  val <- readline(sprintf("  %s (comma-separated, default: %s): ", label, paste(default_vec, collapse = ", ")))
  if (nchar(trimws(val)) == 0) return(default_vec)
  nums <- suppressWarnings(as.numeric(unlist(strsplit(val, "[,;\\s]+"))))
  nums <- nums[!is.na(nums)]
  if (length(nums) != length(default_vec)) {
    cat(sprintf("  Expected %d numbers. Using default.\n", length(default_vec)))
    return(default_vec)
  }
  return(nums)
}

cat(" MODEL PARAMETERS:\n")

params <- list()
params$mcmc      <- input_scalar("mcmc", 10000)
params$burnin    <- input_scalar("burnin", 1000)
params$thinning  <- input_scalar("thinning", 10)
params$py1       <- input_scalar("py1", 0)
params$py2       <- input_scalar("py2", 1.5)
params$pc1       <- input_scalar("pc1", 0)
params$pc2       <- input_scalar("pc2", 2)
params$t1        <- input_scalar("t1", 2)
params$t2        <- input_scalar("t2", min(ntips - 1, 100))
params$scale     <- input_vector("scale", c(20, 10, 5))
params$start     <- input_vector("start", c(1, 1, floor((params$t1 + min(ntips-1, 50))/2)))

params$t2 <- max(params$t1 + 2, min(params$t2, ntips - 1))
params$start[3] <- max(params$t1 + 1, min(params$start[3], params$t2 - 1))

cat(sprintf(" Parameters fixed: t ∈ [%d, %d] | start[3] = %d\n\n", params$t1, params$t2, params$start[3]))

# 7. INTERACTIVE SINGLEPHY TEST LOOP
repeat {
  cat(" Running bgmyc.singlephy on consensus tree...\n")
  
  res_single <- bgmyc.singlephy(
    phylo = tree,
    mcmc = params$mcmc, burnin = params$burnin, thinning = params$thinning,
    py1 = params$py1, py2 = params$py2, pc1 = params$pc1, pc2 = params$pc2,
    t1 = params$t1, t2 = params$t2, scale = params$scale, start = params$start
  )
  
  plot(res_single)
  cat(" Close the graphics window to proceed with convergence diagnostics...\n")
  Sys.sleep(1); flush.console()
  
  #  DYNAMIC CONVERGENCE DIAGNOSTICS 
  ar <- res_single$accept
  cat(sprintf("\n Acceptance rates: py=%.3f | pc=%.3f | th=%.3f\n", ar[1], ar[2], ar[3]))
  
  if (requireNamespace("mcmcse", quietly = TRUE)) {
    ess_vals <- sapply(1:4, function(col) round(mcmcse::ess(res_single$par[, col])))
    cat(sprintf("🔢 ESS: py=%d | pc=%d | th=%d | logL=%d [Optimum: >200]\n",
                ess_vals[1], ess_vals[2], ess_vals[3], ess_vals[4]))
  } else {
    cat("  ESS: package 'mcmcse' not installed. Run install.packages('mcmcse')\n")
    ess_vals <- NULL
  }
  
  recs <- character(0)
  param_names <- c("py", "pc", "th")
  
  if (!is.null(ess_vals) && any(ess_vals < 200)) {
    recs <- c(recs, "• ↑ mcmc or ↓ thinning (ESS < 200: chain needs more independent steps)")
  }
  
  if (length(recs) > 0) {
    cat("\n RECOMMENDATIONS (target range 0.20–0.40):\n")
    for (r in recs) cat(sprintf("   %s\n", r))
    cat("    When re‑running the test, enter suggested values in the scale field.\n")
  } else {
    cat(" Chain convergence stable. Parameters are optimal.\n")
  }

  next_prompt <- if (analysis_mode == 1) {
    "  Proceed to final output (y) or change parameters (n)? [y/n]: "
  } else {
    "  Proceed to multiphylo (y) or change parameters (n)? [y/n]: "
  }
  
  choice <- readline(prompt = sprintf("\n%s", next_prompt))
  if (tolower(choice) != "n") break
  
  cat("\n Update parameters (Enter = keep current):\n")
  
  params$mcmc      <- input_scalar("mcmc", params$mcmc)
  params$burnin    <- input_scalar("burnin", params$burnin)
  params$thinning  <- input_scalar("thinning", params$thinning)
  params$py1       <- input_scalar("py1", params$py1)
  params$py2       <- input_scalar("py2", params$py2)
  params$pc1       <- input_scalar("pc1", params$pc1)
  params$pc2       <- input_scalar("pc2", params$pc2)
  params$t1        <- input_scalar("t1", params$t1)
  params$t2        <- input_scalar("t2", params$t2)
  params$scale     <- input_vector("scale", params$scale)
  params$start     <- input_vector("start", params$start)
  
  if (length(params$scale) != 3) params$scale <- c(20, 10, 5)
  if (length(params$start) != 3) params$start <- c(1, 0.5, floor((params$t1 + params$t2)/2))
  
  params$t2 <- max(params$t1 + 2, min(params$t2, ntips - 1))
  params$start[3] <- max(params$t1 + 1, min(params$start[3], params$t2 - 1))
  
  cat(sprintf(" Parameters updated: t ∈ [%d, %d] | start[3] = %d\n\n", params$t1, params$t2, params$start[3]))
}
# 8. FINAL RESULT PREPARATION
if (analysis_mode == 1) {
  cat("\n Single‑tree analysis finished. Generating outputs...\n")
  final_res <- list(res_single)
  class(final_res) <- "multibgmyc"
} else {
  cat("\n Running bgmyc.multiphylo on selected trees...\n")
  
  n_physical <- parallel::detectCores(logical = FALSE)
  n_workers <- min(n_physical - 1, n_sample)
  if (n_workers < 1) n_workers <- 1
  
  cat(sprintf("   Workers: %d (physical cores: %d)\n", n_workers, n_physical))
  
  plan(multisession, workers = n_workers)
  
  final_res <- future_lapply(seq_along(trees_sample), function(i) {
    bgmyc.singlephy(
      phylo = trees_sample[[i]],
      mcmc = params$mcmc, burnin = params$burnin, thinning = params$thinning,
      py1 = params$py1, py2 = params$py2, pc1 = params$pc1, pc2 = params$pc2,
      t1 = params$t1, t2 = params$t2, scale = params$scale, start = params$start
    )
  }, future.seed = TRUE)
  
  class(final_res) <- "multibgmyc"
  cat(" All trees processed successfully.\n")
  
  #  Gelman‑Rubin (R̂) diagnostics across tree chains
  if (requireNamespace("mcmcse", quietly = TRUE) && length(final_res) > 1) {
    cat("\n Calculating Gelman‑Rubin (R̂) across tree chains...\n")
    
    chains_list <- lapply(final_res, function(res) res$par[, 3])
    names(chains_list) <- paste0("Tree_", seq_along(final_res))
    
    gr_result <- tryCatch(mcmcse::gelman(chains_list), error = function(e) NULL)
    
    if (!is.null(gr_result)) {
      if (!is.null(gr_result$psrf)) {
        rhat <- round(as.numeric(gr_result$psrf[1, "Point est."]), 3)
        cat(sprintf("🔍 Gelman‑Rubin R̂ (threshold): %.3f [Optimum: < 1.05]\n", rhat))
        
        if (rhat > 1.05) {
          cat("  R̂ > 1.05: chains show divergence. Increase mcmc/burnin or check tree topologies.\n")
        } else {
          cat(" Chains converged stably across posterior trees.\n")
        }
      } else {
        cat("  Could not extract R̂ from gelman() result\n")
        cat("   Structure of gr_result:\n")
        print(names(gr_result))
      }
    } else {
      cat("  Could not compute Gelman‑Rubin (possibly too few chains or too short)\n")
    }
  }
}

# 9. VISUALIZATION AND AUTOMATIC CSV EXPORT 
cat("\n Building interactive probability map (CoMa‑style)...\n")
probmat <- spec.probmat(final_res)

# Get tip order from tree
p_tree <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tips_data <- p_tree$data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data$label
n_tips <- length(tip_order)

# Reorder matrix to match tree order
if (!all(rownames(probmat) == tip_order) || !all(colnames(probmat) == tip_order)) {
  probmat <- probmat[tip_order, tip_order]
}

# CUSTOM VISUALIZATION (tree + matrix) 

cat("\n Creating custom visualisation with tree (1:1)...\n")

# If tree_beast not previously loaded, read it via treeio
if (!exists("tree_beast") || !inherits(tree_beast, "treedata")) {
  tree_beast <- treeio::read.beast(consensus_path)
}

# Get tip order from tree
if(!is.null(tree$edge.length)) tree$edge.length[tree$edge.length < 0] <- 0
options(ignore.negative.edge = TRUE)
p_tree_raw <- suppressWarnings(ggtree(tree, layout = "rectangular"))
tree_data_raw <- p_tree_raw$data
tips_data_raw <- tree_data_raw %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data_raw$label
n_tips_tree <- length(tip_order)

cat(sprintf("   Taxa in tree: %d\n", n_tips_tree))
cat(sprintf("   Taxa in probmat: %d\n", nrow(probmat)))

# Find intersection of taxa between tree and matrix
common_taxa <- intersect(tip_order, rownames(probmat))
cat(sprintf("   Common taxa: %d\n", length(common_taxa)))

if (length(common_taxa) == 0) {
  stop(" No common taxa between tree and probability matrix!")
}

# Filter tree to common taxa, preserving annotations via treeio
tips_to_drop <- setdiff(tree_beast@phylo$tip.label, common_taxa)
if(length(tips_to_drop) > 0) {
  tree_beast_filtered <- treeio::drop.tip(tree_beast, tips_to_drop)
} else {
  tree_beast_filtered <- tree_beast
}
tree_filtered <- tree_beast_filtered@phylo

p_tree <- suppressWarnings(ggtree(tree_beast_filtered, layout = "rectangular"))
tree_data <- p_tree$data
tips_data <- tree_data %>% dplyr::filter(isTip) %>% dplyr::arrange(y)
tip_order <- tips_data$label
n_tips <- length(tip_order)

# Filter and reorder matrix
probmat <- probmat[tip_order, tip_order]

cat(sprintf("    Matrix %d×%d synchronised with tree\n", n_tips, n_tips))

# Extract posterior from tree_data annotations
post_col <- intersect(c("posterior", "prob", "Posterior", "PROB"), colnames(tree_data))
if(length(post_col) > 0) {
  post_col <- post_col[1]
  node_posterior <- tree_data[[post_col]]
  names(node_posterior) <- tree_data$node
} else {
  # Fallback: parse from label (if read.nexus left them there)
  node_posterior <- sapply(tree_data$label, function(lbl) {
    if(is.na(lbl) || lbl == "") return(NA_real_)
    m <- regmatches(lbl, regexpr("(posterior|prob)\\s*=\\s*([0-9\\.eE\\-]+)", lbl, ignore.case = TRUE, perl = TRUE))
    if(length(m) == 0 || m == "") return(NA_real_)
    num_str <- sub("^(posterior|prob)\\s*=\\s*", "", m, ignore.case = TRUE, perl = TRUE)
    val <- as.numeric(num_str)
    if(!is.na(val) && val > 1) val <- val / 100
    return(val)
  })
  names(node_posterior) <- tree_data$node
}


# TREE (left panel)

edges <- tree_data %>% dplyr::filter(!is.na(parent))

fig_tree <- plotly::plot_ly()

# Colour function: smooth gradient through several control points
# Red (0) -> Orange (0.25) -> Yellow (0.5) -> Light green (0.75) -> Green (1.0)
get_pp_color <- function(pp) {
  if(is.na(pp)) return("#CCCCCC") # Grey for NA
  pp <- max(0, min(1, pp))
  
  # Define gradient control points
  colors <- list(
    c(0.0, 1.0, 0.0, 0.0),    # Red
    c(0.25, 1.0, 0.5, 0.0),   # Orange
    c(0.5, 1.0, 1.0, 0.0),    # Yellow
    c(0.75, 0.5, 1.0, 0.0),   # Light green
    c(1.0, 0.0, 0.7, 0.0)     # Green
  )
  
  # Find the two nearest control points
  if(pp <= 0) return(sprintf("#%02X%02X%02X", 255, 0, 0))
  if(pp >= 1) return(sprintf("#%02X%02X%02X", 0, 179, 0))
  
  # Interpolate between control points
  for(i in 1:(length(colors)-1)) {
    if(pp >= colors[[i]][1] && pp <= colors[[i+1]][1]) {
      # Linear interpolation between two points
      t <- (pp - colors[[i]][1]) / (colors[[i+1]][1] - colors[[i]][1])
      r <- colors[[i]][2] + t * (colors[[i+1]][2] - colors[[i]][2])
      g <- colors[[i]][3] + t * (colors[[i+1]][3] - colors[[i]][3])
      b <- colors[[i]][4] + t * (colors[[i+1]][4] - colors[[i]][4])
      return(sprintf("#%02X%02X%02X", round(r*255), round(g*255), round(b*255)))
    }
  }
  
  return("#CCCCCC")
}

for(i in 1:nrow(edges)) {
  child <- edges[i, ]
  parent <- tree_data %>% dplyr::filter(node == child$parent)
  if(nrow(parent) == 0) next
  
  # Determine branch colour by posterior of child node
  child_node <- as.character(child$node)
  pp <- node_posterior[child_node]
  branch_color <- get_pp_color(pp)
  
  # Hover text
  pp_text <- ifelse(is.na(pp), "N/A", sprintf("%.2f", pp))
  hover_txt <- paste0("<b>Node:</b> ", child$node, "<br>",
                      "<b>Posterior:</b> ", pp_text)
  
  # Branches (width 5)
  fig_tree <- fig_tree %>% plotly::add_segments(
    x = parent$x, xend = child$x, y = child$y, yend = child$y,
    line = list(color = branch_color, width = 5), hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE)
  fig_tree <- fig_tree %>% plotly::add_segments(
    x = parent$x, xend = parent$x, y = parent$y, yend = child$y,
    line = list(color = branch_color, width = 5), hovertext = hover_txt, hoverinfo = "text", showlegend = FALSE)
}

# Tip labels
max_x <- max(tree_data$x, na.rm = TRUE)
label_offset <- max_x * 0.05 

# Sparse labels 
label_step <- max(1, ceiling(n_tips / 40))
tips_labeled <- tips_data %>%
  dplyr::mutate(row_num = dplyr::row_number()) %>%
  dplyr::filter(row_num %% label_step == 1)

fig_tree <- fig_tree %>% plotly::add_text(
  data = tips_labeled,
  x = ~x + label_offset, y = ~y, text = ~label,
  textposition = "middle left",
  textfont = list(size = 16, family = "monospace", color = "#333333"),
  showlegend = FALSE, hoverinfo = "skip")

# Tree axes layout (autorange = "reversed" so y=1 is on top, like ggtree)
fig_tree <- fig_tree %>% plotly::layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(-0.02, max_x + label_offset + max_x*0.1)),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(0.5, n_tips + 0.5), autorange = "reversed")
)


# HEATMAP (right panel)

hover_text <- matrix(
  paste0("<b>Taxon 1:</b> ", rownames(probmat)[row(probmat)],
         "<br><b>Taxon 2:</b> ", colnames(probmat)[col(probmat)],
         "<br><b>PP(conspecific):</b> ", sprintf("%.2f", probmat)),
  nrow = n_tips)

custom_colorscale <- list(
  list(0.00, "#F0FFFF"),  # Almost white
  list(0.25, "#BDECB6"),  # Light green
  list(0.50, "#98FB98"),  # Medium green
  list(0.95, "#34C924"),  # Dark green
  list(1.00, "#0A5F38")   # Very dark green (maximum confidence)
)

fig_heat <- plotly::plot_ly(
  z = probmat,
  x = 1:n_tips,
  y = 1:n_tips,
  type = "heatmap",
  colorscale = custom_colorscale,  # Or use built‑in viridis palette
  zmin = 0, zmax = 1,
  text = hover_text, hoverinfo = "text",
  showscale = TRUE,
  colorbar = list(
    title = "PP", 
    len = 0.5, 
    x = 1.02,
    tickvals = c(0, 0.25, 0.5, 0.75, 1),
    ticktext = c("0%", "25%", "50%", "75%", "100%")
  )
)

# Heatmap axes (match tree)
fig_heat <- fig_heat %>% plotly::layout(
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(0.5, n_tips + 0.5)),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
               range = c(0.5, n_tips + 0.5), autorange = "reversed")
)


# COMBINE VIA SUBPLOT (1:1)

fig_combined <- plotly::subplot(
  fig_tree, fig_heat,
  nrows = 1,
  widths = c(0.5, 0.5), # Strictly 1:1
  shareY = TRUE,         # Strict Y‑axis synchronisation
  titleX = FALSE, titleY = FALSE
)

fig_combined <- fig_combined %>% plotly::layout(
  title = list(text = "bGMYC4 Conspecificity Probabilities + Consensus Tree", x = 0.5),
  plot_bgcolor = "white",
  hovermode = "closest",
  margin = list(l = 10, r = 60, t = 80, b = 10)
)

# Save
html_path <- file.path(DELIM_DIR, "bGMYC_interactive_heatmap.html")
htmlwidgets::saveWidget(fig_combined, html_path, selfcontained = FALSE, title = "bGMYC Interactive Heatmap")
cat(sprintf(" Saved: %s\n", normalizePath(html_path)))
cat(" Open HTML in browser. Tree on the left, matrix on the right (1:1), all Y‑aligned.\n\n")


# Function to convert cluster list to data frame
create_motu_df <- function(out_list, method_name) {
  motu_list <- lapply(seq_along(out_list), function(i) {
    data.frame(
      Sequence = out_list[[i]],
      MOTU = paste0(method_name, "_", i),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, motu_list)
}

# Automatic delimitation tables
cat("\n Saving delimitation tables...\n")

# Standard threshold (p = 0.05)
out_005 <- bgmyc.point(probmat, ppcutoff = 0.05)
df_005 <- data.frame(
  Sequence = unlist(out_005),
  MOTU_bGMYC = rep(seq_along(out_005), lengths(out_005)),
  stringsAsFactors = FALSE
)
write.table(df_005, file = file.path(DELIM_DIR, "Delimitation_bGMYC_005.csv"), 
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("    p=0.05: Delimitation_bGMYC_005.csv (%d clusters)\n", length(out_005)))

# Strict threshold (p = 0.01)
out_001 <- bgmyc.point(probmat, ppcutoff = 0.01)
df_001 <- data.frame(
  Sequence = unlist(out_001),
  MOTU_bGMYC = rep(seq_along(out_001), lengths(out_001)),
  stringsAsFactors = FALSE
)
write.table(df_001, file = file.path(DELIM_DIR, "Delimitation_bGMYC_001.csv"), 
            row.names = FALSE, sep = ";", dec = ".", quote = FALSE, fileEncoding = "UTF-8")
cat(sprintf("    p=0.01: Delimitation_bGMYC_001.csv (%d clusters)\n", length(out_001)))

# Export full probability table
spec_out  <- bgmyc.spec(final_res)
write.csv(spec_out$specprobs, file.path(DELIM_DIR, "bGMYC_delimitation_results.csv"), 
          row.names = FALSE, fileEncoding = "UTF-8")
cat("    Probabilities: bGMYC_delimitation_results.csv\n")

cat("\n Top 10 clusters (p=0.05):\n")
for (i in seq_along(out_005)[1:min(10, length(out_005))]) {
  taxa  <- out_005[[i]]
  cat(sprintf("   %2d: %s%s\n", i, paste(taxa[1:min(3, length(taxa))], collapse = ", "),
              if (length(taxa) > 3) "..." else ""))
}

cat("\n Interactive bGMYC4 analysis complete!\n")