---
layout: default
title: "Bosmina Integrative Delimitation"
description: "Supplementary materials for integrative species delimitation of Bosmina"
---

# 🧬 Bosmina Integrative Delimitation
**Supplementary Materials & Interactive Visualizations**

Welcome to the supplementary materials repository for integrative species delimitation analysis of *Bosmina* (Cladocera: Bosminidae). This repository contains all data, scripts, and results supporting our research.

---

## 📂 Repository Structure

<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 25px; margin: 40px 0;">

<div style="border: 3px solid #1F78B4; border-radius: 15px; padding: 25px; background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%); box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
<h3 style="color: #1F78B4; margin-top: 0; display: flex; align-items: center;">📊 Raw Data</h3>
<p><strong>Location:</strong> <code>/data/</code></p>
<p>Original delimitation results, sequence alignments, and input files for all analyses.</p>
<ul style="margin-bottom: 20px;">
<li>Delimitation results per method (CSV)</li>
<li>Sequence alignments (FASTA)</li>
<li>Phylogenetic trees (Newick)</li>
<li>Metadata and annotations</li>
</ul>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/data" style="display: inline-block; background: #1F78B4; color: white; padding: 12px 24px; text-decoration: none; border-radius: 8px; font-weight: bold; text-align: center; width: 100%; box-sizing: border-box;">Browse Raw Data →</a>
</div>

<div style="border: 3px solid #E41A1C; border-radius: 15px; padding: 25px; background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%); box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
<h3 style="color: #E41A1C; margin-top: 0; display: flex; align-items: center;">💻 Scripts</h3>
<p><strong>Location:</strong> <code>/scripts/</code></p>
<p>Complete source code for all computational analyses and visualization pipelines.</p>
<ul style="margin-bottom: 20px;">
<li>R scripts for data processing</li>
<li>Agreement matrix calculation</li>
<li>Phylogenetic visualization</li>
<li>Statistical analysis pipelines</li>
</ul>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/scripts" style="display: inline-block; background: #E41A1C; color: white; padding: 12px 24px; text-decoration: none; border-radius: 8px; font-weight: bold; text-align: center; width: 100%; box-sizing: border-box;">View Scripts →</a>
</div>

<div style="border: 3px solid #4DAF4A; border-radius: 15px; padding: 25px; background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%); box-shadow: 0 4px 6px rgba(0,0,0,0.1);">
<h3 style="color: #4DAF4A; margin-top: 0; display: flex; align-items: center;">📈 Results</h3>
<p><strong>Location:</strong> <code>/results/</code></p>
<p>Comprehensive analysis outputs including tables, figures, and interactive visualizations.</p>
<ul style="margin-bottom: 20px;">
<li>Delimitation summaries</li>
<li>Statistical reports</li>
<li>Phylogenetic trees</li>
<li><strong>Interactive visualizations</strong></li>
</ul>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/results" style="display: inline-block; background: #4DAF4A; color: white; padding: 12px 24px; text-decoration: none; border-radius: 8px; font-weight: bold; text-align: center; width: 100%; box-sizing: border-box;">Explore Results →</a>
</div>

</div>

---

## 🌟 Featured Interactive Visualization

<div style="border: 4px solid #984EA3; border-radius: 20px; padding: 30px; background: linear-gradient(135deg, #f5f3ff 0%, #ede9fe 100%); box-shadow: 0 8px 16px rgba(152,78,163,0.2); margin: 40px 0;">

<h2 style="color: #984EA3; margin-top: 0;">🎨 Agreement Matrix + Clade-colored Phylogeny</h2>

<p style="font-size: 1.1em; line-height: 1.6;">
Explore the congruence between multiple species delimitation methods through our interactive heatmap visualization. The matrix is synchronized with a phylogenetic tree colored by major clades.
</p>

<div style="background: white; padding: 20px; border-radius: 10px; margin: 20px 0;">
<h4 style="margin-top: 0; color: #666;">✨ Features:</h4>
<ul style="line-height: 1.8;">
<li><strong>Hover tooltips</strong> showing taxon names, agreement percentages, and method counts</li>
<li><strong>Color-coded clades</strong>: <span style="color: #E41A1C;">■</span> Eubosmina, <span style="color: #377EB8;">■</span> Liederobosmina, <span style="color: #4DAF4A;">■</span> Lunobosmina, <span style="color: #984EA3;">■</span> Colombian clade</li>
<li><strong>Fully zoomable and interactive</strong> interface</li>
<li><strong>Lightweight HTML</strong> (~5 MB) for fast loading</li>
</ul>
</div>

<a href="results/interactive/Delimitation_heatmap_bgmyc_tree.html" style="display: inline-block; background: #984EA3; color: white; padding: 15px 40px; text-decoration: none; border-radius: 10px; font-weight: bold; font-size: 1.2em; box-shadow: 0 4px 8px rgba(152,78,163,0.3);">🚀 Open Interactive Visualization</a>

</div>

---

## 📊 Results Breakdown

<div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 40px 0;">

<div style="border: 2px solid #FF7F00; border-radius: 12px; padding: 20px; background: #fff;">
<h4 style="color: #FF7F00; margin-top: 0;">📋 Delimitation Results</h4>
<p>Summary tables and congruence matrices for all delimitation methods.</p>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/results/delimitation" style="color: #FF7F00; font-weight: bold; text-decoration: none;">View Results →</a>
</div>

<div style="border: 2px solid #A65628; border-radius: 12px; padding: 20px; background: #fff;">
<h4 style="color: #A65628; margin-top: 0;">📉 Statistics</h4>
<p>Statistical summaries, congruence scores, and method comparisons.</p>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/results/statistics" style="color: #A65628; font-weight: bold; text-decoration: none;">View Statistics →</a>
</div>

<div style="border: 2px solid #F781BF; border-radius: 12px; padding: 20px; background: #fff;">
<h4 style="color: #F781BF; margin-top: 0;">🌳 Phylogenetic Trees</h4>
<p>Newick and Nexus tree files with clade annotations.</p>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/results/trees" style="color: #F781BF; font-weight: bold; text-decoration: none;">View Trees →</a>
</div>

<div style="border: 2px solid #999999; border-radius: 12px; padding: 20px; background: #fff;">
<h4 style="color: #666; margin-top: 0;">🎭 Interactive Figures</h4>
<p>Plotly-based interactive visualizations and dashboards.</p>
<a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/tree/main/results/interactive" style="color: #666; font-weight: bold; text-decoration: none;">Browse Interactive →</a>
</div>

</div>

---

## 📚 About This Study

| **Organisms** | *Bosmina* spp. (Cladocera: Bosminidae) |
|---------------|----------------------------------------|
| **Focus** | Integrative species delimitation |
| **Clades** | *Eubosmina*, *Liederobosmina*, *Lunobosmina*, Colombian clade |
| **Methods** | ABGD, ASAP, bPTP, GMYC, PTP, and others |
| **License** | [GNU GPL v3.0](https://github.com/DmitryKarabanov/Bosmina_delimitaions/blob/main/LICENSE) |

---

## 👨‍🔬 Contact & Citation

**Dmitry Karabanov**

For questions about the data, scripts, or analysis, please [open an issue on GitHub](https://github.com/DmitryKarabanov/Bosmina_delimitaions/issues).

### How to Cite

If you use these data or scripts in your research, please cite:

> Karabanov D. et al. (2026). Integrative species delimitation of *Bosmina*... *[Journal name, in press]*

---

<p style="text-align: center; color: #666; font-size: 0.9em; margin-top: 50px; padding: 20px; border-top: 2px solid #eee;">
<em>This supplementary material is openly available under the <a href="https://github.com/DmitryKarabanov/Bosmina_delimitaions/blob/main/LICENSE">GNU GPL v3.0 License</a>. Last updated: 2026.</em>
</p>
