#!/usr/bin/env python3
"""
Publication-ready phylogenetic tree visualization for Hypsizygus marmoreus
mitochondrial 15-gene concatenated tree, colored by geographic origin.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib import font_manager
from Bio import Phylo
from io import StringIO
import numpy as np
import re

# --- CJK font setup ---
_cjk_font = '/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc'
import os
if os.path.exists(_cjk_font):
    font_manager.fontManager.addfont(_cjk_font)
    prop = font_manager.FontProperties(fname=_cjk_font)
    plt.rcParams['font.sans-serif'] = [prop.get_name()] + plt.rcParams.get('font.sans-serif', [])
    plt.rcParams['axes.unicode_minus'] = False
    print(f"Using CJK font: {prop.get_name()}")
else:
    print("Warning: CJK font not found")

# ============================================================
# 1. Read tree
# ============================================================
tree_file = "04_phylogeny/mito_15genes_partition.treefile"
tree = Phylo.read(tree_file, "newick")

# ============================================================
# 2. Sample metadata: geographic origin
# ============================================================
location_map = {
    # GenBank reference genomes
    "MF133443.1":  "Fuzhou, China",
    "MH382825.1":  "Korea",
    # SRA - Korea
    "SRR7874787":  "Korea",
    # SRA - Fujian, China (separate study strains)
    "SRR12151860": "Fujian, China",
    "SRR12151871": "Fujian, China",
    "SRR12151875": "Fujian, China",
    "SRR12151883": "Fujian, China",
    # SRA - Fuzhou, China (FQX_MS01 population)
    "SRR8699796":  "Fuzhou, China",
    "SRR8699797":  "Fuzhou, China",
    "SRR8699800":  "Fuzhou, China",
    "SRR8699801":  "Fuzhou, China",
    "SRR8699802":  "Fuzhou, China",
    "SRR8699803":  "Fuzhou, China",
    "SRR8699804":  "Fuzhou, China",
    "SRR8699805":  "Fuzhou, China",
    "SRR8699808":  "Fuzhou, China",
    "SRR8699809":  "Fuzhou, China",
    "SRR8699811":  "Fuzhou, China",
    "SRR8699813":  "Fuzhou, China",
    "SRR8699814":  "Fuzhou, China",
    "SRR8699815":  "Fuzhou, China",
    "SRR8699816":  "Fuzhou, China",
    "SRR8699817":  "Fuzhou, China",
    "SRR8699833":  "Fuzhou, China",
    "SRR8699834":  "Fuzhou, China",
    "SRR8699835":  "Fuzhou, China",
    "SRR8699837":  "Fuzhou, China",
    # User's own samples - Shanghai
    "f2":          "Shanghai, China",
    "f4":          "Shanghai, China",
    "nn12-1":      "Shanghai, China",
    "nn12-17":     "Shanghai, China",
}

# Strain / display name mapping
display_name_map = {
    "MF133443.1":  "MF133443.1",
    "MH382825.1":  "MH382825.1",
    "SRR7874787":  "SRR7874787 (51987-8)",
    "SRR12151860": "SRR12151860 (Hyz5-9)",
    "SRR12151871": "SRR12151871 (Hyz5-78)",
    "SRR12151875": "SRR12151875 (HM-Hy46)",
    "SRR12151883": "SRR12151883 (Hm-G2)",
    "f2":          "f2",
    "f4":          "f4",
    "nn12-1":      "nn12-1",
    "nn12-17":     "nn12-17",
}

# Color palette for locations
color_map = {
    "Shanghai, China": "#E63946",   # red
    "Fujian, China":   "#457B9D",   # steel blue
    "Fuzhou, China":   "#2A9D8F",   # teal
    "Korea":           "#E9C46A",   # gold
}

# ============================================================
# 3. Helper: extract bootstrap from confidence attribute
# ============================================================
def get_bootstrap(clade):
    """Extract bootstrap value from clade confidence or name."""
    if clade.confidence is not None:
        return clade.confidence
    if clade.name and '/' in str(clade.name):
        # format: "SH-aLRT/UFBoot" e.g. "73.4/65"
        parts = str(clade.name).split('/')
        try:
            return float(parts[-1])
        except:
            pass
    if clade.name:
        try:
            return float(clade.name)
        except:
            pass
    return None

# ============================================================
# 4. Tree layout calculation (rectangular cladogram/phylogram)
# ============================================================
def get_all_terminals(clade):
    """Recursively get all terminal names."""
    if clade.is_terminal():
        return [clade.name]
    names = []
    for child in clade.clades:
        names.extend(get_all_terminals(child))
    return names

def assign_y_positions(clade, y_counter, y_positions):
    """Assign y positions to terminals first, then internal nodes."""
    if clade.is_terminal():
        y_positions[id(clade)] = y_counter[0]
        y_counter[0] += 1
    else:
        for child in clade.clades:
            assign_y_positions(child, y_counter, y_positions)
        child_ys = [y_positions[id(c)] for c in clade.clades]
        y_positions[id(clade)] = np.mean(child_ys)

def assign_x_positions(clade, x_start, x_positions, use_branch_length=True):
    """Assign x positions based on branch lengths."""
    x_positions[id(clade)] = x_start
    for child in clade.clades:
        bl = child.branch_length if (use_branch_length and child.branch_length) else 0.001
        assign_x_positions(child, x_start + bl, x_positions, use_branch_length)

# ============================================================
# 5. Draw the tree
# ============================================================
fig, ax = plt.subplots(1, 1, figsize=(14, 13))
fig.patch.set_facecolor('white')

# Layout
y_positions = {}
x_positions = {}
assign_y_positions(tree.root, [0], y_positions)
assign_x_positions(tree.root, 0, x_positions, use_branch_length=True)

# Get max x for scaling
max_x = max(x_positions.values())
if max_x == 0:
    max_x = 1

# Track which vertical connectors already drawn to avoid redundancy
_drawn_verticals = set()

def draw_tree(clade, ax):
    """Recursively draw the tree."""
    x = x_positions[id(clade)]
    y = y_positions[id(clade)]

    if not clade.is_terminal() and len(clade.clades) >= 2:
        # Draw vertical connector ONCE per internal node
        child_ys = [y_positions[id(c)] for c in clade.clades]
        y_min, y_max = min(child_ys), max(child_ys)
        if id(clade) not in _drawn_verticals:
            ax.plot([x, x], [y_min, y_max],
                    color='#333333', linewidth=1.0, solid_capstyle='round')
            _drawn_verticals.add(id(clade))

    for child in clade.clades:
        cx = x_positions[id(child)]
        cy = y_positions[id(child)]

        # Horizontal line (branch)
        ax.plot([x, cx], [cy, cy], color='#333333', linewidth=1.0, solid_capstyle='round')

        # Bootstrap label on internal branches
        if not child.is_terminal():
            bs = get_bootstrap(child)
            if bs is not None and bs >= 65:
                ax.text(cx + max_x * 0.002, cy + 0.35, f'{int(bs)}',
                        fontsize=6.5, color='#C0392B', ha='left', va='bottom',
                        fontweight='bold')

        draw_tree(child, ax)

    # Terminal label
    if clade.is_terminal():
        name = clade.name if clade.name else ""
        loc = location_map.get(name, "Unknown")
        color = color_map.get(loc, "#888888")
        display = display_name_map.get(name, name)

        # Dashed line from tip to label alignment position
        label_x = max_x * 1.02
        ax.plot([x, label_x], [y, y], color='#dddddd', linewidth=0.4,
                linestyle=':', zorder=1)

        # Colored dot
        ax.scatter(label_x + max_x * 0.005, y, s=50, color=color, zorder=5,
                   edgecolors='#444444', linewidths=0.4)
        # Label text
        ax.text(label_x + max_x * 0.018, y, display,
                fontsize=8, va='center', ha='left', color='#222222',
                fontfamily='sans-serif')

draw_tree(tree.root, ax)

# ============================================================
# 6. Legend & styling
# ============================================================
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map["Shanghai, China"],
           markeredgecolor='#444444', markeredgewidth=0.4,
           markersize=10, label='Shanghai, China (上海)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map["Fujian, China"],
           markeredgecolor='#444444', markeredgewidth=0.4,
           markersize=10, label='Fujian, China (福建)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map["Fuzhou, China"],
           markeredgecolor='#444444', markeredgewidth=0.4,
           markersize=10, label='Fuzhou, China (福州)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color_map["Korea"],
           markeredgecolor='#444444', markeredgewidth=0.4,
           markersize=10, label='Korea (韩国)'),
]

leg = ax.legend(handles=legend_elements, loc='upper left', frameon=True,
                fontsize=9, framealpha=0.95, edgecolor='#bbbbbb',
                title='Geographic Origin', title_fontsize=10,
                fancybox=True, shadow=False,
                bbox_to_anchor=(0.0, 1.0))
leg.get_title().set_fontweight('bold')

# Scale bar
scale_len = 0.005  # substitutions per site
n_tips = len(list(tree.get_terminals()))
scale_y = -2.5
ax.plot([0, scale_len], [scale_y, scale_y], color='black', linewidth=2,
        solid_capstyle='butt')
# Ticks at ends
for sx in [0, scale_len]:
    ax.plot([sx, sx], [scale_y - 0.3, scale_y + 0.3], color='black', linewidth=1.5)
ax.text(scale_len / 2, scale_y - 1.0, f'{scale_len} subs/site', fontsize=8,
        ha='center', va='top', fontstyle='italic')

# Title
ax.set_title('Phylogenetic tree of $\\it{H. marmoreus}$ mitochondrial genomes\n'
             '(15 protein-coding genes, partition model, IQ-TREE 3.0)',
             fontsize=13, fontweight='bold', pad=18)

# Clean up axes
ax.set_xlim(-max_x * 0.02, max_x * 1.55)
ax.set_ylim(-5, n_tips + 1)
ax.axis('off')

plt.tight_layout()

# Save
output_png = "04_phylogeny/mito_phylogeny_publication.png"
output_pdf = "04_phylogeny/mito_phylogeny_publication.pdf"
plt.savefig(output_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(output_pdf, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_png}")
print(f"Saved: {output_pdf}")
plt.close()
