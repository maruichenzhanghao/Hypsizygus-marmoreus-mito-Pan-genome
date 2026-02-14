# Hypsizygus marmoreus Mitochondrial Pan-genome

## 真姬菇（*Hypsizygus marmoreus*）线粒体泛基因组分析

[![GitHub](https://img.shields.io/badge/GitHub-Hypsizygus--marmoreus--mito--Pan--genome-blue)](https://github.com/maruichenzhanghao/Hypsizygus-marmoreus-mito-Pan-genome)

---

## 📖 项目概述

本项目对 **31 株真姬菇** 线粒体基因组进行了泛基因组分析，包括：
- 线粒体基因组组装与注释
- 基于 Cactus 的泛基因组图构建
- 变异检测与统计分析
- 系统发育分析
- dN/dS 选择压力分析
- 内含子多态性分析
- 覆盖深度验证
- 基因重排分析

## 📊 主要发现

| 分析模块 | 关键结果 |
|----------|----------|
| **泛基因组图** | 217 节点 / 293 边 / 总长 220,364 bp |
| **变异统计** | 2,506 个变异: SNPs=1,563, InDels=522, MNPs=246 |
| **区域分布** | 基因间区 78.9%, CDS 15.6%, rRNA 5.4%, tRNA ~0% |
| **选择压力** | 15 个核心基因全部 ω < 1 → 纯化选择 |
| **内含子** | cox1 有 5 种内含子模式, 内含子多态性丰富 |
| **基因重排** | IR 介导的 2 个倒位事件 |

## 🗂️ 仓库结构

```
.
├── 01.assembly/                    # 线粒体基因组组装
│   ├── batch_mito_analysis.sh      # 批量组装流程
│   └── stats_plots/                # 组装统计与可视化
│
├── 02.annotation/                  # 基因组注释
│   ├── gb/                         # GenBank 格式注释脚本
│   │   ├── extract_cds_for_pan_genome.py
│   │   ├── pangenome.quxian.py     # 泛基因组曲线
│   │   ├── pangenome.bingtu.py     # 饼图可视化
│   │   ├── pangenome.yanbenzhuzhuangtu.py  # 样本柱状图
│   │   └── cunzaiqueshiretu.py     # 基因存在/缺失热图
│   └── trna/                       # tRNA 注释
│
├── 03.pangenome/                   # 泛基因组构建 (Cactus)
│   ├── mito_samples.txt            # 样本列表
│   ├── generate_input.sh           # 生成 Cactus 输入
│   ├── run_cactus.sh               # Cactus 运行脚本
│   ├── run_cactus1.sh              # Cactus 步骤1
│   ├── run_cactus2.sh              # Cactus 步骤2
│   ├── run_cactus_3.sh             # Cactus 步骤3
│   └── run_odgi.sh                 # ODGI 可视化
│
├── 04.variant_analysis/            # 变异检测与分析
│   ├── vcf_cell/                   # VCF 变异统计
│   │   ├── mito_gb_to_gff.py       # GB → GFF 转换
│   │   ├── calculate_variant_density.py
│   │   ├── gene_variant_stats.py
│   │   └── publication_ready_charts.py
│   └── pangenome_results/          # 泛基因组分析结果
│       ├── analyze_graph_pangenome.py
│       ├── corrected_graph_analysis.py
│       ├── deep_graph_analysis.py
│       ├── plot_pangenome_variants.py
│       ├── plot_combined_figure.py
│       └── plot_horizontal.py
│
├── 05.phylogeny/                   # 系统发育分析
│   ├── build_phylogeny.py          # 建树流程
│   ├── plot_tree.py                # 树的可视化
│   └── data/                       # 基因序列与比对
│
├── 06.dnds/                        # dN/dS 选择压力
│   ├── plot_dnds_ggplot.R          # R 可视化
│   ├── dnds_data_for_R.csv         # 数据
│   └── dnds_summary.tsv            # 汇总表
│
├── 07.intron_analysis/             # 内含子分析
│   ├── plotD_intron_heatmap.R
│   ├── plotE_gene_structure.R
│   ├── plotF_intron_genome.R
│   ├── plotG_intron_orf_pangenome.R
│   └── data/                       # 内含子统计数据
│
├── 08.coverage_analysis/           # 覆盖深度分析
│   ├── run_coverage.sh
│   └── analyze_coverage.py
│
├── 09.rearrangement/               # 基因重排分析
│   ├── plot_IR_mechanism.py
│   ├── plot_IR_mechanism_v2.py
│   ├── plot_rearrangement.R
│   └── tongji_sv.sh
│
├── figures/                        # 发表用图片
├── README.md
└── .gitignore
```

## 🔧 依赖软件

### 生物信息学工具
| 软件 | 版本 | 用途 |
|------|------|------|
| Cactus | ≥ 2.6 | 泛基因组图构建 |
| ODGI | - | 图基因组可视化 |
| vg | - | 变异图操作 |
| IQ-TREE | 2 | 系统发育建树 |
| MAFFT | ≥ 7 | 序列比对 |
| SAMtools | ≥ 1.15 | BAM 文件处理 |
| Minimap2 | ≥ 2.24 | 序列比对 |
| BCFtools | ≥ 1.15 | VCF 处理 |

### Python 环境
```bash
pip install biopython matplotlib seaborn pandas numpy scipy networkx gffutils
```

### R 包
```R
install.packages(c("ggplot2", "ggtree", "pheatmap", "dplyr", "tidyr", "reshape2"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
```

## 🚀 分析流程

### 1. 线粒体基因组组装
```bash
# 使用 MitoHiFi / GetOrganelle 等工具从 HiFi/Illumina 数据组装
bash 01.assembly/batch_mito_analysis.sh
```

### 2. 基因组注释
```bash
# 使用 MITOS / MFAnnot 进行自动注释，人工校正
# GenBank 格式注释文件存放在 02.annotation/gb/
```

### 3. Cactus 泛基因组构建
```bash
# 生成输入文件
bash 03.pangenome/generate_input.sh
# 运行 Cactus 多步流程
bash 03.pangenome/run_cactus.sh
bash 03.pangenome/run_cactus1.sh
bash 03.pangenome/run_cactus2.sh
bash 03.pangenome/run_cactus_3.sh
# ODGI 可视化
bash 03.pangenome/run_odgi.sh
```

### 4. 变异检测与统计
```bash
python 04.variant_analysis/vcf_cell/gene_variant_stats.py
python 04.variant_analysis/vcf_cell/calculate_variant_density.py
python 04.variant_analysis/vcf_cell/publication_ready_charts.py
```

### 5. 系统发育分析
```bash
python 05.phylogeny/build_phylogeny.py   # 15 基因串联 + IQ-TREE
python 05.phylogeny/plot_tree.py         # 可视化
```

### 6. dN/dS 分析
```bash
Rscript 06.dnds/plot_dnds_ggplot.R
```

### 7. 内含子分析
```bash
Rscript 07.intron_analysis/plotD_intron_heatmap.R
Rscript 07.intron_analysis/plotE_gene_structure.R
Rscript 07.intron_analysis/plotF_intron_genome.R
Rscript 07.intron_analysis/plotG_intron_orf_pangenome.R
```

### 8. 覆盖深度验证
```bash
bash 08.coverage_analysis/run_coverage.sh
python 08.coverage_analysis/analyze_coverage.py
```

## 📋 样本信息

本研究使用 **31 株真姬菇** 线粒体基因组，包括：
- **4 株** 自测 HiFi 序列 (f2, f4, nn12-1, nn12-17)
- **2 株** NCBI 参考基因组 (MF133443.1, MH382825.1)
- **25 株** 公共数据库 SRA 数据 (SRR 系列)

## 📄 引用

如果本项目对你有帮助，请引用：

> [论文信息待补充]

## 📝 License

本项目采用 MIT License 开源协议。

## 📬 联系方式

如有任何问题，请通过 [GitHub Issues](https://github.com/maruichenzhanghao/Hypsizygus-marmoreus-mito-Pan-genome/issues) 联系。
