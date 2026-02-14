# Mitochondrial Graph-Based Pan-Genome Analysis of *Hypsizygus marmoreus*

**Structural Variation, Adaptive Evolution and Its Implications for Germplasm Resource Improvement**

Ruichen Ma<sup>1†</sup>, Wenyun Li<sup>2†</sup>, Yongmei Miao<sup>1</sup>, Ruiheng Yang<sup>2,3</sup>, Youran Shao<sup>2</sup>, Junjun Shang<sup>2,3</sup>, Yan Li<sup>2,3</sup>, Yuan Gao<sup>2</sup>, Dapeng Bao<sup>2,3\*</sup>, Yingying Wu<sup>2,3\*</sup>

<sup>1</sup> Anhui Science and Technology University, Fengyang 233100, China  
<sup>2</sup> National Engineering Research Center of Edible Fungi, Institute of Edible Fungi, Shanghai Academy of Agricultural Sciences, Shanghai 201403, China  
<sup>3</sup> Shanghai Key Laboratory of Agricultural Genetics and Breeding, Shanghai Academy of Agricultural Sciences, Shanghai 201106, China  

<sup>†</sup> These authors contributed equally to this work.  
<sup>\*</sup> Correspondence: baodapeng@saas.sh.cn (D.B.); wuyingying@sibs.ac.cn (Y.W.)

---

## Abstract

As semi-autonomous organelles, mitochondria function through the coordinated regulation of nuclear genomes and their own genetic material. In this study, we conducted a comparative mitochondrial genome analysis of **31 *Hypsizygus marmoreus* strains** (4 newly sequenced monokaryons and 27 public datasets). The mitochondrial genome sizes ranged from 98,284 to 111,087 bp, exhibiting significant structural diversity driven by non-coding region dynamics and intronic polymorphisms. The 31 mitochondrial genomes were assembled into a **graph-based pan-genome (220,364 bp, 217 nodes)** capturing abundant SNPs, InDels, and structural variations. Eight gene rearrangement patterns and five genetic clusters were identified, providing breeding-relevant genetic markers and a genomic framework for germplasm classification and molecular breeding of *H. marmoreus*.

**Keywords:** *Hypsizygus marmoreus*; Mitochondria; Structural variation; Graph-based pan-genome; Germplasm resource development

---

## Figure Materials

Each figure directory contains the plotting scripts used to generate the corresponding figure, along with a README describing the methods, tools, and panel descriptions.

| Figure | Description | Directory |
|--------|-------------|-----------|
| **Fig. 1** | Mitochondrial genome map (circular, OGDRAW) | [`Figure1/`](Figure1/) |
| **Fig. 2** | tRNA secondary structure prediction | [`Figure2/`](Figure2/) |
| **Fig. 3** | Repetitive sequence heatmap (SSR, interspersed, tandem) | [`Figure3/`](Figure3/) |
| **Fig. 4** | Relative synonymous codon usage (RSCU) analysis | [`Figure4/`](Figure4/) |
| **Fig. 5** | Ka/Ks selective pressure analysis (15 core PCGs) | [`Figure5/`](Figure5/) |
| **Fig. 6** | Nucleotide diversity (π) and gene rearrangement patterns | [`Figure6/`](Figure6/) |
| **Fig. 7** | Phylogenetic tree and population structure (ADMIXTURE + IQ-TREE) | [`Figure7/`](Figure7/) |
| **Fig. 8** | Pan-genome gene families, intron analysis, and PAV | [`Figure8/`](Figure8/) |
| **Fig. 9** | Graph-based pan-genome construction and structural variation | [`Figure9/`](Figure9/) |

### Figure Overview

<p align="center">
  <img src="figures/Figure1.png" width="45%" alt="Figure 1"/>
  <img src="figures/Figure7.png" width="45%" alt="Figure 7"/>
</p>
<p align="center">
  <img src="figures/Figure8.png" width="45%" alt="Figure 8"/>
  <img src="figures/Figure9.png" width="45%" alt="Figure 9"/>
</p>

---

## Workflows and Methods

### Genome Assembly

Scripts for mitochondrial genome assembly from HiFi and Illumina data.  
→ [`workflows/Assembly/`](workflows/Assembly/)

| Tool | Version | Purpose |
|------|---------|---------|
| Flye | v2.9.5 | *De novo* HiFi assembly |
| GetOrganelle | v1.7.7.1 | Illumina mitochondrial assembly |
| Bandage | v0.9.0 | Assembly graph visualization |

### Genome Annotation

Scripts for mitochondrial genome annotation and format conversion.  
→ [`workflows/Annotation/`](workflows/Annotation/)

| Tool | Version | Purpose |
|------|---------|---------|
| MFannot | — | Automated mitochondrial annotation |
| Mitos | — | Automated mitochondrial annotation |
| Geneious | v2025.0.2 | Manual annotation correction |
| OGDRAW | v1.3.1 | Circular genome map visualization |
| tRNAscan-SE | v2.0.12 | tRNA gene prediction |

### Graph Pan-Genome Construction

Scripts for Cactus-based graph pan-genome construction and variant calling.  
→ [`workflows/Pangenome_Construction/`](workflows/Pangenome_Construction/)

| Tool | Version | Purpose |
|------|---------|---------|
| Minigraph-Cactus | v2.9.3 | Graph pan-genome construction |
| vg toolkit | v1.6.1.0 | Variant calling (`vg call`) |
| bcftools | v1.19 | VCF statistics and filtering |
| ODGI | v0.9.2 | Graph topology analysis |
| Bandage | v0.9.0 | GFA graph visualization |

### Comparative Analysis

| Tool | Version | Purpose |
|------|---------|---------|
| MAFFT | v7.526 | Multiple sequence alignment |
| IQ-TREE | v2.0.7 | Maximum-likelihood phylogeny |
| ADMIXTURE | v1.3.0 | Population structure analysis |
| pamlX | v1.3.1 | Ka/Ks calculation |
| DnaSP | v6 | SNP detection and π calculation |
| OrthoFinder | v2.55 | Gene family clustering |
| PhyloSuite | v1.2.3 | PCG sequence extraction |
| MEGA | v11 | RSCU calculation |
| MISA | — | SSR identification |
| REPuter | — | Interspersed repeat detection |
| TRF | — | Tandem repeat detection |

---

## Sample Information

This study used **31 *H. marmoreus* strains**:

| Category | Strains | Source |
|----------|---------|--------|
| Self-sequenced HiFi monokaryons | f2, f4, nn12-1, nn12-17 | This study (PV946885, PX600725, PX600726, PX600727) |
| NCBI reference genomes | MF133443.1, MH382825.1 | GenBank |
| Public WGS datasets | 25 strains (SRR series) | NCBI SRA |

---

## Key Findings

| Analysis | Result |
|----------|--------|
| Genome size range | 98,284 – 111,087 bp |
| Graph pan-genome | 220,364 bp, 217 nodes, 293 edges |
| Total variants | 2,506 sites: SNPs (1,563), InDels (522), MNPs (246), Others (175) |
| Variant distribution | Intergenic 78.9%, CDS 15.6%, rRNA 5.4%, tRNA ~0% |
| Selective pressure | All 15 core genes: ω < 1 (purifying selection) |
| Rearrangement | 8 gene arrangement patterns, 2 IR-mediated inversion events |
| Population structure | 5 genetic clusters (K=5 optimal) |

---

## Repository Structure

```
.
├── README.md
├── figures/                           # Publication figures (PNG)
│   ├── Figure1.png ~ Figure9.png
│
├── Figure1/                           # Mitochondrial genome map
├── Figure2/                           # tRNA secondary structure
├── Figure3/                           # Repetitive sequence heatmap
├── Figure4/                           # RSCU codon analysis
│   └── mimazi.test1.R
├── Figure5/                           # Ka/Ks selective pressure
│   ├── kaks_15gene.R
│   ├── KAKS_xiangxiantu.R
│   └── plot_dnds_ggplot.R
├── Figure6/                           # Nucleotide diversity + gene order
│   ├── gene_pi_plot_english.R
│   └── geneorder.duose.R
├── Figure7/                           # Phylogeny + population structure
│   ├── build_phylogeny.py
│   ├── plot_tree.py
│   ├── admixture.R / admixture.2.R / admixture.3.R
│   ├── CVerror.R
│   └── PCA3D.R
├── Figure8/                           # Pan-genome + intron analysis
│   ├── pangenome.quxian.py            # Pan-genome curve
│   ├── pangenome.bingtu.py            # Pie chart
│   ├── cunzaiqueshiretu.py            # PAV heatmap
│   ├── COX1.R                         # Gene structure
│   └── plotD~G_intron_*.R             # Intron analysis plots
├── Figure9/                           # Graph pan-genome + SV
│   ├── plot_pangenome_variants.py
│   ├── plot_IR_mechanism_v2.py
│   ├── plot_rearrangement.R
│   ├── ggtree_sv_tree.R
│   └── analyze_coverage.py
│
└── workflows/                         # Analysis pipelines
    ├── Assembly/                      # Genome assembly scripts
    ├── Annotation/                    # Genome annotation scripts
    └── Pangenome_Construction/        # Cactus pan-genome pipeline
```

## Environment Setup

### Python
```bash
pip install biopython matplotlib seaborn pandas numpy scipy networkx
```

### R
```R
install.packages(c("ggplot2", "pheatmap", "ggsci", "ggh4x", "plotly",
                    "reshape2", "dplyr", "tidyr", "ape", "rgl"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio", "gggenes"))
```

---

## Citation

> Ma R, Li W, Miao Y, Yang R, Shao Y, Shang J, Li Y, Gao Y, Bao D\*, Wu Y\*. Mitochondrial Graph-Based Pan-Genome Analysis of *Hypsizygus marmoreus*: Structural Variation, Adaptive Evolution and Its Implications for Germplasm Resource Improvement. *[Journal]*, 2025.

---

## License

This project is licensed under the MIT License.

## Contact

For questions, please open a [GitHub Issue](https://github.com/maruichenzhanghao/Hypsizygus-marmoreus-mito-Pan-genome/issues).
