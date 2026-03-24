"""
Elbow plot + PC count comparison for the dimensionality reduction blog series.
Blog post: quantavo.ai/blog/dimensionality-reduction-part-2
Author: Cory Henn | quantavo.ai

Uses the 10x Genomics PBMC 3k dataset (downloaded automatically by scanpy).
Produces a 3-panel figure: elbow plot, UMAP at 15 PCs, UMAP at 100 PCs.
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

adata = sc.datasets.pbmc3k_processed()

# Recompute PCA with enough components to show the full elbow
sc.tl.pca(adata, n_comps=100, random_state=42)
var_ratio = adata.uns['pca']['variance_ratio']

# UMAP from 15 PCs (good choice — at the elbow)
sc.pp.neighbors(adata, n_pcs=15, random_state=42, key_added='good')
sc.tl.umap(adata, random_state=42, neighbors_key='good')
umap_good = adata.obsm['X_umap'].copy()

# UMAP from 100 PCs (too many — noise leaks in)
sc.pp.neighbors(adata, n_pcs=100, random_state=42, key_added='bad')
sc.tl.umap(adata, random_state=42, neighbors_key='bad')
umap_bad = adata.obsm['X_umap'].copy()

# Plot
BG = '#1c2b38'
FG = '#dde5e8'
GRD = '#3a5060'

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(28, 8), facecolor=BG)
fig.subplots_adjust(wspace=0.3)
fig.suptitle('How many PCs should you keep?',
             fontsize=18, color=FG, fontweight='bold', y=1.02)

# Panel 1: elbow plot
ax1.set_facecolor(BG)
ax1.plot(range(1, 51), var_ratio[:50], '-o', color='#06b6d4',
         markersize=4, linewidth=2, alpha=0.9)
ax1.axvline(x=15, color='#f97316', linestyle='--', linewidth=2, alpha=0.8)
ax1.text(16.5, var_ratio[2], '15 PCs', color='#f97316',
         fontsize=13, fontweight='bold')
ax1.fill_betweenx([0, var_ratio[0]], 0, 15, alpha=0.08, color='#06b6d4')
ax1.text(7, var_ratio[5], 'signal', color='#06b6d4',
         fontsize=12, ha='center', style='italic', alpha=0.7)
ax1.text(35, var_ratio[20], 'noise', color='#ef4444',
         fontsize=12, ha='center', style='italic', alpha=0.7)
ax1.set_title('Elbow Plot', fontsize=16, fontweight='bold', color=FG, pad=12)
ax1.set_xlabel('Principal Component', color=FG, fontsize=11)
ax1.set_ylabel('Variance Ratio', color=FG, fontsize=11)
ax1.tick_params(colors=GRD, labelsize=9)
for spine in ax1.spines.values():
    spine.set_edgecolor(GRD)

# Panels 2 and 3: UMAP comparisons
umap_panels = [
    (ax2, umap_good, '15 PCs \u2014 clean separation',
     'Signal preserved', '#55A868'),
    (ax3, umap_bad, '100 PCs \u2014 noise corrupts the structure',
     'Noise leaks in', '#C44E52'),
]

for ax, coords, title, note, note_color in umap_panels:
    ax.set_facecolor(BG)
    for ct in adata.obs['louvain'].cat.categories:
        mask = adata.obs['louvain'] == ct
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   s=10, alpha=0.6, linewidths=0, label=ct)
    ax.set_title(title, fontsize=16, fontweight='bold', color=FG, pad=12)
    ax.set_xlabel('UMAP 1', color=FG, fontsize=11)
    ax.set_ylabel('UMAP 2', color=FG, fontsize=11)
    ax.tick_params(colors=GRD, labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRD)
    ax.text(0.97, 0.04, note, transform=ax.transAxes, ha='right', va='bottom',
            fontsize=11, color=note_color, style='italic', fontweight='semibold')
    legend = ax.legend(fontsize=7, frameon=True, loc='upper right',
                       framealpha=0.3, edgecolor=GRD, markerscale=2, ncol=2)
    legend.get_frame().set_facecolor(BG)
    for text in legend.get_texts():
        text.set_color(FG)

fig.text(0.98, 0.003, 'quantavo.ai  |  Data: 10x Genomics PBMC 3k',
         ha='right', fontsize=9, color='#666', style='italic')

plt.tight_layout()
fig.savefig('elbow-pc-comparison.png', dpi=200, bbox_inches='tight',
            facecolor=BG, edgecolor='none')
plt.close()
print('Saved: elbow-pc-comparison.png')
