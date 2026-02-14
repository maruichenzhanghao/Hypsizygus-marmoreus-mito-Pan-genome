#!/usr/bin/env python3
"""
IR-mediated mitochondrial genome rearrangement mechanism — V2 compact figure
Based on actual f2 reference genome coordinates (106,605 bp, GenBank).

Panel (a): Event A — 385 bp IR → 18.5 kb inversion (nad2–nad3–rps3)
Panel (b): Event B — DR → IR → 12.2 kb inversion (atp8–nad5–nad4L–cox3)

Actual coordinates from f2.gb / f2.mito.clean.gff / f2_self_repeats.txt:
  Event A IR: f2:15,994–16,373 ↔ f2:34,145–34,529 (RC), 385 bp, 88.1%
  Event B DR: f2:41,639–41,902 ↔ f2:77,358–77,615, 265 bp, 90.6%
  nad2: 18,789–20,654 (+)
  nad3: 20,654–21,040 (+)
  rps3: 27,694–29,484 (-)
  atp8: 67,498–67,656 (+)
  nad5: 71,614–74,971 (-)
  nad4L: 74,971–75,237 (-)
  cox3: 76,569–77,378 (-)
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8,
    'axes.linewidth': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'savefig.dpi': 500,
})

# ── Colors ──
C_IR     = '#C0392B'
C_IR2    = '#E74C3C'
C_DR     = '#8E44AD'
C_FWD    = '#2471A3'    # blue: forward strand genes
C_REV    = '#27AE60'    # green: inverted genes
C_TRNA   = '#F5B041'    # tRNA clusters
C_FLANK  = '#ABB2B9'
C_RECOM  = '#E67E22'
C_BG     = '#F8F9F9'

# ── Helper functions ──
def gene_arrow(ax, x, y, w, h, label, color, direction='+', fs=6.5, ec='#2C3E50'):
    """Gene arrow polygon."""
    tip = 0.15 * w
    if direction == '+':
        pts = [(x, y-h/2), (x+w-tip, y-h/2), (x+w, y), (x+w-tip, y+h/2), (x, y+h/2)]
    else:
        pts = [(x+w, y-h/2), (x+tip, y-h/2), (x, y), (x+tip, y+h/2), (x+w, y+h/2)]
    poly = plt.Polygon(pts, closed=True, fc=color, ec=ec, lw=0.5, zorder=3, alpha=0.92)
    ax.add_patch(poly)
    ax.text(x+w/2, y, label, ha='center', va='center', fontsize=fs,
            fontweight='bold', fontstyle='italic', color='white', zorder=4,
            path_effects=[pe.withStroke(linewidth=1.2, foreground='#2C3E50')])

def ir_block(ax, x, y, w, h, label, color=C_IR, fs=5.5):
    """IR/DR element."""
    rect = FancyBboxPatch((x, y-h/2), w, h, boxstyle="round,pad=0.015",
                           fc=color, ec='#1B2631', lw=0.6, zorder=3, alpha=0.95)
    ax.add_patch(rect)
    ax.text(x+w/2, y, label, ha='center', va='center', fontsize=fs,
            fontweight='bold', color='white', zorder=4)

def trna_block(ax, x, y, w, h, label='tRNAs', color=C_TRNA, fs=5):
    """tRNA cluster block."""
    rect = FancyBboxPatch((x, y-h/2), w, h, boxstyle="round,pad=0.01",
                           fc=color, ec='#7D6608', lw=0.4, zorder=3, alpha=0.7)
    ax.add_patch(rect)
    ax.text(x+w/2, y, label, ha='center', va='center', fontsize=fs,
            color='#7D6608', fontstyle='italic', zorder=4)

def flank_line(ax, x1, x2, y, color=C_FLANK, lw=3):
    """Flanking region as thick line."""
    ax.plot([x1, x2], [y, y], color=color, lw=lw, solid_capstyle='round', zorder=1)

def coord_label(ax, x, y, text, fs=5, color='#566573'):
    """Position coordinate label."""
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, color=color, fontstyle='italic')


# ============================================================
# FIGURE — compact: ~170 mm × 100 mm
# ============================================================
fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(6.7, 4.2),
                                  gridspec_kw={'height_ratios': [1, 1.15], 'hspace': 0.35})

for ax in (ax_a, ax_b):
    ax.set_xlim(-0.5, 18.5)
    ax.set_ylim(-3.8, 3.8)
    ax.set_aspect('auto')
    ax.axis('off')

gh = 0.55  # gene height

# ============================================================
# Panel (a) — Event A
# ============================================================
ax_a.text(-0.3, 3.6, '(a)', fontsize=11, fontweight='bold', va='top')
ax_a.text(0.5, 3.6, 'Event A: 385 bp IR-mediated inversion (18,535 bp)',
          fontsize=8.5, fontweight='bold', va='top', color='#1B2631')

# ── Row 1: Reference (f2) ──
y1 = 2.0
ax_a.text(-0.3, y1, 'Ref.', fontsize=7, ha='center', va='center',
          fontweight='bold', color='#2C3E50')

flank_line(ax_a, 0.3, 1.2, y1)
ir_block(ax_a, 1.3, y1, 0.85, gh, 'IR-A', C_IR)
# trnG cluster
trna_block(ax_a, 2.3, y1, 0.7, gh*0.8, 'trnG')
# nad2(+) — f2:18,789–20,654
gene_arrow(ax_a, 3.15, y1, 1.9, gh, 'nad2', C_FWD, '+')
# nad3(+) — f2:20,654–21,040
gene_arrow(ax_a, 5.2, y1, 0.8, gh, 'nad3', C_FWD, '+', fs=5.5)
# tRNA cluster (trnM, trnS, trnP, trnN, trnL, trnE, trnS, trnH)
trna_block(ax_a, 6.15, y1, 1.6, gh*0.8, '8 tRNAs')
# rps3(-) — f2:27,694–29,484
gene_arrow(ax_a, 7.9, y1, 1.9, gh, 'rps3', C_FWD, '-')
# orf region
trna_block(ax_a, 10.0, y1, 0.7, gh*0.8, 'orfs', '#D5DBDB', fs=5)
# IR-A' (RC)
ir_block(ax_a, 10.85, y1, 0.85, gh, "IR-A'", C_IR2)
flank_line(ax_a, 11.85, 12.7, y1)

# Coordinate labels
coord_label(ax_a, 1.72, y1 - 0.53, '15,994', fs=4.5)
coord_label(ax_a, 11.27, y1 - 0.53, '34,529', fs=4.5)

# 18,535 bp bracket
ax_a.annotate('', xy=(10.7, y1 + 0.6), xytext=(2.15, y1 + 0.6),
              arrowprops=dict(arrowstyle='<->', color='#2C3E50', lw=0.8))
ax_a.text(6.4, y1 + 0.8, '18,535 bp', fontsize=6, ha='center', color='#2C3E50')

# IR identity box
ax_a.text(14.2, y1, '385 bp\n88.1% identity', fontsize=6, ha='center', va='center',
          color=C_IR, fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.2', fc='#FDEDEC', ec=C_IR, lw=0.6, alpha=0.9))

# ── Recombination arrow ──
y_mid = 0.0
# Cross lines representing recombination
ax_a.annotate('', xy=(10.5, y_mid+0.15), xytext=(2.2, y1-0.45),
              arrowprops=dict(arrowstyle='->', color=C_RECOM, lw=1.2,
                              connectionstyle='arc3,rad=0.2'))
ax_a.annotate('', xy=(2.2, y_mid+0.15), xytext=(10.5, y1-0.45),
              arrowprops=dict(arrowstyle='->', color=C_RECOM, lw=1.2,
                              connectionstyle='arc3,rad=-0.2'))
ax_a.text(13, y_mid+0.15, 'Intramolecular\nhomologous\nrecombination',
          fontsize=5.5, ha='center', va='center', color=C_RECOM,
          fontweight='bold', fontstyle='italic')

# ── Row 2: Inverted ──
y2 = -2.0
ax_a.text(-0.3, y2, 'Inv.', fontsize=7, ha='center', va='center',
          fontweight='bold', color='#E74C3C')

flank_line(ax_a, 0.3, 1.2, y2)
ir_block(ax_a, 1.3, y2, 0.85, gh, 'IR-A', C_IR)
trna_block(ax_a, 2.3, y2, 0.7, gh*0.8, 'orfs', '#D5DBDB', fs=5)
# INVERTED: rps3(+), tRNAs reversed, nad3(-), nad2(-)
gene_arrow(ax_a, 3.15, y2, 1.9, gh, 'rps3', C_REV, '+')
trna_block(ax_a, 5.2, y2, 1.6, gh*0.8, '8 tRNAs')
gene_arrow(ax_a, 6.95, y2, 0.8, gh, 'nad3', C_REV, '-', fs=5.5)
gene_arrow(ax_a, 7.9, y2, 1.9, gh, 'nad2', C_REV, '-')
trna_block(ax_a, 9.95, y2, 0.7, gh*0.8, 'trnG')
ir_block(ax_a, 10.85, y2, 0.85, gh, "IR-A'", C_IR2)
flank_line(ax_a, 11.85, 12.7, y2)

# Inversion bar
ax_a.annotate('', xy=(10.7, y2 - 0.6), xytext=(2.15, y2 - 0.6),
              arrowprops=dict(arrowstyle='<->', color='#E74C3C', lw=0.8))
ax_a.text(6.4, y2 - 0.8, '180° inversion', fontsize=6, ha='center',
          color='#E74C3C', fontstyle='italic')

# Sample count box
ax_a.text(14.2, y2, '12/31 strains\n(38.7%)', fontsize=6, ha='center', va='center',
          color='#2C3E50', fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.2', fc='#EBF5FB', ec='#2980B9', lw=0.6))

# Thin separator
ax_a.axhline(y=-3.6, xmin=0.02, xmax=0.98, color='#D5D8DC', lw=0.4, ls='--')


# ============================================================
# Panel (b) — Event B: DR → IR → Inversion
# ============================================================
ax_b.text(-0.3, 3.6, '(b)', fontsize=11, fontweight='bold', va='top')
ax_b.text(0.5, 3.6, 'Event B: Dispersed repeat → IR → inversion (12,171 bp)',
          fontsize=8.5, fontweight='bold', va='top', color='#1B2631')

# ── Row 1: f2 reference with DR (same direction) ──
y1b = 2.2
ax_b.text(-0.3, y1b, 'f2', fontsize=7, ha='center', va='center',
          fontweight='bold', color=C_DR)

flank_line(ax_b, 0.3, 0.9, y1b)
gene_arrow(ax_b, 1.0, y1b, 0.9, gh, 'atp6', C_FWD, '-', fs=5.5)
# DR1 in dpo region
ir_block(ax_b, 2.05, y1b, 0.8, gh, 'DR₁', C_DR, fs=5)
gene_arrow(ax_b, 3.0, y1b, 0.65, gh, 'dpo', '#7F8C8D', '-', fs=5)

# rnl region
gene_arrow(ax_b, 3.8, y1b, 1.4, gh*0.85, 'rnl', '#7F8C8D', '-', fs=6)

# trnR(-), trnC(-)
trna_block(ax_b, 5.35, y1b, 0.55, gh*0.7, 'R,C', C_TRNA, fs=4)
# atp8(+)
gene_arrow(ax_b, 6.05, y1b, 0.65, gh, 'atp8', C_FWD, '+', fs=5)
# tRNAs (W, L, R)
trna_block(ax_b, 6.85, y1b, 0.55, gh*0.7, 'W,L,R', C_TRNA, fs=4)
# nad5(-)
gene_arrow(ax_b, 7.55, y1b, 1.5, gh, 'nad5', C_FWD, '-', fs=6)
# nad4L(-)
gene_arrow(ax_b, 9.2, y1b, 0.7, gh, 'nad4L', C_FWD, '-', fs=5)
# cox3(-)
gene_arrow(ax_b, 10.05, y1b, 0.9, gh, 'cox3', C_FWD, '-', fs=5.5)
# DR2 near cox3-orf389
ir_block(ax_b, 11.1, y1b, 0.8, gh, 'DR₂', C_DR, fs=5)
# orf389
trna_block(ax_b, 12.05, y1b, 0.55, gh*0.7, 'orfs', '#D5DBDB', fs=4)
flank_line(ax_b, 12.75, 13.4, y1b)

# DR arrows (same direction →)
ax_b.annotate('→', xy=(2.45, y1b+0.48), fontsize=7, ha='center', color=C_DR, fontweight='bold')
ax_b.annotate('→', xy=(11.5, y1b+0.48), fontsize=7, ha='center', color=C_DR, fontweight='bold')
ax_b.text(6.8, y1b+0.85, 'Direct repeats (same orientation)', fontsize=5, ha='center',
          color=C_DR, fontstyle='italic')

# DR info
ax_b.text(15.0, y1b, '265 bp DR\n90.6% identity', fontsize=5.5, ha='center', va='center',
          color=C_DR, fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.2', fc='#F5EEF8', ec=C_DR, lw=0.6))

# ── Step arrow ──
y_step = 0.85
ax_b.annotate('', xy=(6.5, y_step-0.25), xytext=(6.5, y_step+0.35),
              arrowprops=dict(arrowstyle='->', color='#566573', lw=1.5))
ax_b.text(9.5, y_step+0.05, 'One copy inverts → forms IR pair',
          fontsize=5.5, ha='center', va='center', color='#566573', fontstyle='italic')

# ── Row 2: After DR→IR conversion ──
y2b = 0.0
ax_b.text(-0.3, y2b, 'IR\nformed', fontsize=5.5, ha='center', va='center',
          fontweight='bold', color=C_IR)

flank_line(ax_b, 0.3, 0.9, y2b)
gene_arrow(ax_b, 1.0, y2b, 0.9, gh, 'atp6', C_FWD, '-', fs=5.5)
ir_block(ax_b, 2.05, y2b, 0.8, gh, 'IR-B', C_IR, fs=5)
gene_arrow(ax_b, 3.0, y2b, 1.4, gh*0.85, 'rnl', '#7F8C8D', '-', fs=6)
trna_block(ax_b, 4.55, y2b, 0.55, gh*0.7, 'R,C', C_TRNA, fs=4)
gene_arrow(ax_b, 5.25, y2b, 0.65, gh, 'atp8', C_FWD, '+', fs=5)
trna_block(ax_b, 6.05, y2b, 0.55, gh*0.7, 'W,L,R', C_TRNA, fs=4)
gene_arrow(ax_b, 6.75, y2b, 1.5, gh, 'nad5', C_FWD, '-', fs=6)
gene_arrow(ax_b, 8.4, y2b, 0.7, gh, 'nad4L', C_FWD, '-', fs=5)
gene_arrow(ax_b, 9.25, y2b, 0.9, gh, 'cox3', C_FWD, '-', fs=5.5)
ir_block(ax_b, 10.3, y2b, 0.8, gh, "IR-B'", C_IR2, fs=5)
trna_block(ax_b, 11.25, y2b, 0.55, gh*0.7, 'orfs', '#D5DBDB', fs=4)
flank_line(ax_b, 11.95, 12.6, y2b)

# IR arrows (opposite direction)
ax_b.annotate('→', xy=(2.45, y2b+0.48), fontsize=7, ha='center', color=C_IR, fontweight='bold')
ax_b.annotate('←', xy=(10.7, y2b+0.48), fontsize=7, ha='center', color=C_IR2, fontweight='bold')

# Recombination cross
ax_b.annotate('', xy=(10.0, y2b-0.5-0.9), xytext=(2.85, y2b-0.4),
              arrowprops=dict(arrowstyle='->', color=C_RECOM, lw=1.0,
                              connectionstyle='arc3,rad=0.2'))
ax_b.annotate('', xy=(2.85, y2b-0.5-0.9), xytext=(10.0, y2b-0.4),
              arrowprops=dict(arrowstyle='->', color=C_RECOM, lw=1.0,
                              connectionstyle='arc3,rad=-0.2'))
ax_b.text(13.5, y2b - 0.7, 'IR-mediated\nrecombination',
          fontsize=5.5, ha='center', va='center', color=C_RECOM,
          fontweight='bold', fontstyle='italic')

# IR info
ax_b.text(15.0, y2b, '213 bp IR\n92.0% identity', fontsize=5.5, ha='center', va='center',
          color=C_IR, fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.2', fc='#FDEDEC', ec=C_IR, lw=0.6))

# ── Row 3: Inverted product ──
y3b = -2.2
ax_b.text(-0.3, y3b, 'Inv.', fontsize=7, ha='center', va='center',
          fontweight='bold', color='#E74C3C')

flank_line(ax_b, 0.3, 0.9, y3b)
gene_arrow(ax_b, 1.0, y3b, 0.9, gh, 'atp6', C_FWD, '-', fs=5.5)
ir_block(ax_b, 2.05, y3b, 0.8, gh, 'IR-B', C_IR, fs=5)
gene_arrow(ax_b, 3.0, y3b, 1.4, gh*0.85, 'rnl', '#7F8C8D', '+', fs=6)
# INVERTED genes: cox3(+), nad4L(+), nad5(+), atp8(-)
gene_arrow(ax_b, 4.55, y3b, 0.9, gh, 'cox3', C_REV, '+', fs=5.5)
gene_arrow(ax_b, 5.6, y3b, 0.7, gh, 'nad4L', C_REV, '+', fs=5)
gene_arrow(ax_b, 6.45, y3b, 1.5, gh, 'nad5', C_REV, '+', fs=6)
trna_block(ax_b, 8.1, y3b, 0.55, gh*0.7, 'R,L,W', C_TRNA, fs=4)
gene_arrow(ax_b, 8.8, y3b, 0.65, gh, 'atp8', C_REV, '-', fs=5)
trna_block(ax_b, 9.6, y3b, 0.55, gh*0.7, 'C,R', C_TRNA, fs=4)
ir_block(ax_b, 10.3, y3b, 0.8, gh, "IR-B'", C_IR2, fs=5)
trna_block(ax_b, 11.25, y3b, 0.55, gh*0.7, 'orfs', '#D5DBDB', fs=4)
flank_line(ax_b, 11.95, 12.6, y3b)

# Inversion bar
ax_b.annotate('', xy=(10.1, y3b - 0.6), xytext=(3.0, y3b - 0.6),
              arrowprops=dict(arrowstyle='<->', color='#E74C3C', lw=0.8))
ax_b.text(6.55, y3b - 0.82, '12,171 bp inversion (180° flip)', fontsize=5.5,
          ha='center', color='#E74C3C', fontstyle='italic')

# Sample count
ax_b.text(15.0, y3b, '5/31 strains\n(16.1%)', fontsize=6, ha='center', va='center',
          color='#2C3E50', fontweight='bold',
          bbox=dict(boxstyle='round,pad=0.2', fc='#EBF5FB', ec='#2980B9', lw=0.6))

# ── Legend (compact, at bottom) ──
leg_y = -3.55
items = [
    (C_IR, 's', 'Inverted repeat (IR)'),
    (C_FWD, '>', 'Reference orientation'),
    (C_REV, '>', 'Inverted orientation'),
    (C_DR, 's', 'Dispersed repeat (DR)'),
    (C_RECOM, 'X', 'Recombination'),
]
for i, (col, mk, lab) in enumerate(items):
    lx = 0.3 + i * 3.5
    ax_b.plot(lx, leg_y, mk, color=col, markersize=5, zorder=5)
    ax_b.text(lx+0.25, leg_y, lab, fontsize=5, va='center', color='#2C3E50')

# ============================================================
# Save
# ============================================================
out_png = '/home/maxinxin/workspace/01.mito/03.pangenome.cactus/mito_results815/dpo_check/IR_rearrangement_mechanism_v2.png'
out_pdf = '/home/maxinxin/workspace/01.mito/03.pangenome.cactus/mito_results815/dpo_check/IR_rearrangement_mechanism_v2.pdf'

fig.savefig(out_png, dpi=500, bbox_inches='tight', facecolor='white', pad_inches=0.15)
fig.savefig(out_pdf, bbox_inches='tight', facecolor='white', pad_inches=0.15)
plt.close()
print(f"✓ V2 saved:\n  PNG: {out_png}\n  PDF: {out_pdf}")
