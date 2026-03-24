"""
Synthetic scRNA-seq data + figure for the dimensionality reduction blog series.
Blog post: quantavo.ai/blog/dimensionality-reduction-part-1
Author: Cory Henn | quantavo.ai

600 cells (5 immune types x 120), 200 genes. Compares raw 2-gene
projection vs PCA on top HVGs.
"""

import numpy as np
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.random.seed(0)

N_CELLS_PER_TYPE = 120
N_GENES = 200
CELL_TYPES = ['T cells', 'B cells', 'Macrophages', 'NK cells', 'Dendritic cells']
N_TYPES = len(CELL_TYPES)
PALETTE = ['#4C72B0', '#DD8452', '#C44E52', '#64B5CD', '#55A868']

# Each type gets 20 marker genes (high expression) in genes 100-199.
# Genes 0-99 are housekeeping — genes 7 and 83 are both in this range.
BASELINE = 2.0
MARKER_EXP = 15.0
NOISE_STD = 2.5

cell_means = []
for i in range(N_TYPES):
    means = np.ones(N_GENES) * BASELINE
    means[100 + i * 20 : 100 + i * 20 + 20] = MARKER_EXP
    cell_means.append(means)

X_list, y_list = [], []
for i, means in enumerate(cell_means):
    cells = means + np.random.randn(N_CELLS_PER_TYPE, N_GENES) * NOISE_STD
    X_list.append(cells)
    y_list.extend([i] * N_CELLS_PER_TYPE)

X = np.vstack(X_list)
y = np.array(y_list)

# Center but don't standardize — keeps natural variance magnitudes
X_centered = X - X.mean(axis=0)

# Top 80 HVGs by variance, then PCA
hvg_idx = np.argsort(X_centered.var(axis=0))[-80:]
pca = PCA(n_components=2, random_state=42)
X_pca = pca.fit_transform(X_centered[:, hvg_idx])

# Raw panel: genes 7 and 83 (both housekeeping, no cell-type signal)
X_raw = X_centered[:, [7, 83]]

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))
fig.subplots_adjust(wspace=0.32)
fig.suptitle('600 cells, 200 genes \u2014 same data, different views',
             color='#555', fontsize=11, y=1.01)

panels = [
    (ax1, X_raw, 'Gene 7 expression', 'Gene 83 expression',
     'Raw High-Dimensional Projection', 'No structure visible', '#C44E52'),
    (ax2, X_pca, 'PC 1', 'PC 2',
     'After Dimensionality Reduction', 'Cell types clearly separate', '#55A868'),
]

for ax, coords, xlab, ylab, title, note, note_color in panels:
    ax.set_xlabel(xlab, fontsize=10)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=10)

    for i, ct in enumerate(CELL_TYPES):
        mask = y == i
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=PALETTE[i], s=22, alpha=0.80, linewidths=0, label=ct)

    ax.text(0.97, 0.04, note, transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color=note_color, style='italic', fontweight='semibold')

ax2.legend(fontsize=9, frameon=True, loc='upper right',
           framealpha=0.8, markerscale=1.5)

fig.savefig('dim-reduction-comparison.png', dpi=180, bbox_inches='tight')
plt.close(fig)

print(f"PC1: {pca.explained_variance_ratio_[0]:.1%}")
print(f"PC2: {pca.explained_variance_ratio_[1]:.1%}")
print("Saved: dim-reduction-comparison.png")
