# Figure 9: Graph-based pan-genome and structural variation analysis

## Description

- **(A)** Construction of the mitochondrial graph pan-genome of *H. marmoreus* (Bandage visualization of GFA)
- **(B)** Variant statistics:
  - (a) Variant type distribution: 1,563 SNPs, 522 InDels, 246 MNPs, 175 others
  - (b) Variant distribution across genomic regions (intergenic 78.9%, CDS 15.6%, rRNA 5.4%)
  - (c) Mutation spectrum (Ts/Tv = 1.19)
- **(C)** IR-mediated gene rearrangement mechanism
- **(D)** SV-based phylogenetic tree (NJ tree with ggtree)
- Coverage depth validation across 28 samples (4 HiFi + 24 Illumina)

## Scripts

| Script | Description |
|--------|-------------|
| `plot_pangenome_variants.py` | Combined variant type/region/spectrum figures |
| `plot_IR_mechanism_v2.py` | IR-mediated rearrangement mechanism diagram |
| `plot_rearrangement.R` | Gene rearrangement visualization (R version) |
| `ggtree_sv_tree.R` | SV-based NJ phylogenetic tree with ggtree |
| `plot_combined_figure.py` | Combined multi-panel figure |
| `plot_horizontal.py` | Horizontal layout variant figure |
| `analyze_coverage.py` | Coverage depth analysis and validation |

## Tools Used

- **Minigraph-Cactus v2.9.3** — Graph pan-genome construction
- **vg toolkit v1.6.1.0** — Variant calling (`vg call`)
- **bcftools v1.19** — VCF statistics
- **Bandage v0.9.0** — GFA visualization
- **ODGI v0.9.2** — Graph topology analysis
- **R (ggtree, ggplot2)** — Tree and plot visualization
- **Python (matplotlib, seaborn, networkx)** — Visualization

![Figure9](../figures/Figure9.png)
