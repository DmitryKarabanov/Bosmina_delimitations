---
layout: default
title: Interactive Visualizations
---


#  Interactive Visualizations

This directory contains all interactive HTML visualizations supporting the integrative species delimitation analysis of *Bosmina* (Cladocera: Bosminidae).
Each visualization is a self-contained HTML file that can be opened directly in any modern web browser. No special software or internet connection is required.

---

##  Available Visualizations

| # | Visualization | Main file | Description |
|---|---------------|------|-------------|
| 1 |  **Net** | `Bosmina_TCS_MedianJoining` | A nice haplotype network for the entire sample of Bosmina  |
| 2 |  **bGMYC4** | `bGMYC_interactive_heatmap` | bGMYC4 results with phylogenetic tree and heatmap |
| 3 |  **Summary** | `Delimitation_heatmap_bgmyc_tree` | Comparison charts of all delimitation methods |
| 4 |  **Agreement** | `Agreement_Matrix` | Interactive heatmap showing congruence between delimitation methods, synchronized with clade-colored phylogenetic tree |


---

##  Featured Visualization

Our main interactive figure combines two synchronized views:

###  Left Panel: Clade-Colored Phylogeny
- **Branches colored** by major clades 
- **Hover** any branch to see taxon name / clade assignment / branch length

###  Right Panel: Agreement Matrix
- **Heatmap** showing pairwise agreement between delimitation method(s)
- **Color scale:** drom 0% agreement to 100% agreement
- **Hover** any cell to see taxon names / agreement percentage 

### 🔹 Synchronization
Both panels share the same Y-axis, so each row in the matrix corresponds exactly to a tip on the tree.

---

##  How to Use

### Option 1: View on GitHub Pages (recommended)
Visit: **[Net](https://dmitrykarabanov.github.io/Bosmina_delimitaions/interactive/Bosmina_TCS_MedianJoining.html)**

Visit: **[bGMYC4](https://dmitrykarabanov.github.io/Bosmina_delimitaions/interactive/bGMYC_interactive_heatmap.html)**

Visit: **[Summary](https://dmitrykarabanov.github.io/Bosmina_delimitaions/interactive/Delimitation_heatmap_bgmyc_tree.html)**

Visit: **[Agreement](https://dmitrykarabanov.github.io/Bosmina_delimitaions/interactive/Agreement_Matrix.html)**


### Option 2: Open locally
1. Download the `.html` file and the additional directory
2. Double-click to open in your web browser
3. Use mouse to **zoom**, **pan**, and **hover** for details

### Option 3: Embed in presentations
- Take screenshots for static figures
- Use browser "Print → Save as PDF" for publication-ready outputs
- Export Plotly data via the camera icon in the toolbar

---

##  Troubleshooting

**Q: The visualization doesn't load.**  
A: Try a different browser or clear your cache. Some corporate networks block large HTML files.

**Q: Hover tooltips don't appear.**  
A: Make sure JavaScript is enabled in your browser settings.

**Q: The file is too large to view on GitHub.**  
A: Download the file locally — GitHub has size limits for inline HTML rendering.

---

##  Related Resources

-  **[Raw Data](../data/)** — Input files used to generate these visualizations
-  **[Scripts](../scripts/)** — Source code (R) for reproducing all figures
-  **[Results](../results/)** — Static tables and statistical outputs


[← Back to Main Page](../README.md)
