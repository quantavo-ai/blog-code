"""
t-SNE vs UMAP comparison for the dimensionality reduction blog series.
Blog post: quantavo.ai/blog/dimensionality-reduction-part-3
Author: Cory Henn | quantavo.ai

Uses the 10x Genomics PBMC 3k dataset (downloaded automatically by scanpy).
Both embeddings computed from the top 30 of 50 PCs for a fair comparison.
"""

import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

adata = sc.datasets.pbmc3k_processed()

# Recompute PCA, neighbors, and both embeddings from the same PCs
sc.tl.pca(adata, n_comps=50, random_state=42)
sc.pp.neighbors(adata, n_pcs=30, random_state=42)
sc.tl.tsne(adata, random_state=42, n_pcs=30, perplexity=30)
sc.tl.umap(adata, random_state=42)

tsne_coords = adata.obsm['X_tsne'].copy()
umap_coords = adata.obsm['X_umap'].copy()

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 18))
fig.suptitle('2,638 PBMCs - same cells, same PCs, different methods',
             fontsize=14, color='#333', fontweight='bold', y=0.99)

panels = [
    (ax1, tsne_coords, 't-SNE', 't-SNE 1', 't-SNE 2', 'upper right'),
    (ax2, umap_coords, 'UMAP', 'UMAP 1', 'UMAP 2', 'upper left'),
]

for ax, coords, title, xlab, ylab, legend_loc in panels:
    for ct in adata.obs['louvain'].cat.categories:
        mask = adata.obs['louvain'] == ct
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   s=14, alpha=0.6, linewidths=0, label=ct)
    ax.set_title(title, fontsize=18, fontweight='bold', pad=12)
    ax.set_xlabel(xlab, fontsize=12)
    ax.set_ylabel(ylab, fontsize=12)
    ax.legend(fontsize=8, framealpha=0.8, loc=legend_loc,
              markerscale=2.5, ncol=2)

fig.text(0.98, 0.003, 'quantavo.ai  |  Data: 10x Genomics PBMC 3k',
         ha='right', fontsize=9, color='#999', style='italic')

plt.tight_layout()
fig.savefig('tsne-vs-umap-comparison.png', dpi=200, bbox_inches='tight')
plt.close()
print('Saved: tsne-vs-umap-comparison.png')
