#  Analysis Scripts

Complete source code for all computational analyses.

## Scripts

### [Delimitation](delimitaion/)
- **`bGMYC_interactive.R`** — Agreement matrix calculation and interactive visualization (Plotly + ggtree)
- **`locMin_delimitation.R`** — Calculation of local minima (locMin), two distance variants
- ...

### [Genetic diversity and Net](genepop/)
- Scripts for calculating genetic diversity and constructing an interactive haplotype network

## Requirements

- R >= 4.5
- Packages: `ape`, `dplyr`, `tidyr`, ets.
- Bioconductor packages for phylogenetic analysis
- Python >= 3.1
- - Packages: `numpy`, `scipy`, `pyvis`, ets.

## Usage
For bGMYC delimitation you need our bGMYC4 package and R>=4.5
```r
bGMYC4_interactive.R
```
For nice haplotype network use Anaconda / Python>=3.1
```Python
network_new.py
```

## Output
Scripts generate results in /results/ directory.

---

##  Related Resources

-  **[Raw Data](../data/)** — Input files used to generate these visualizations
-  **[Interactive](../interactive/)** — Source code (R) for reproducing all figures
-  **[Results](../results/)** — Static tables and statistical outputs

[← Back to Main Page](../README.md)


_This supplementary material is openly available under the GNU GPL v3.0 License. Last updated: 2026._
