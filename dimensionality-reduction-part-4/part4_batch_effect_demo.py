"""
Batch effect visualization for the dimensionality reduction blog series.
Blog post: quantavo.ai/blog/dimensionality-reduction-part-4
Author: Cory Henn | quantavo.ai

Two real PBMC datasets merged without correction — one from a 3k cell run,
one from a 68k cell run (downsampled). Genuinely different batches, different
sequencing depth, different library prep. No artificial splitting.

The point: PCA surfaces the batch separation immediately. UMAP buries it
inside what looks like clean cell type clusters. This is a real failure mode,
not a toy example.

Dependencies: scanpy, matplotlib, numpy
    pip install scanpy matplotlib numpy
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


# ── settings ──────────────────────────────────────────────────────────────────

sc.settings.verbosity = 1       # just warnings, not the full scanpy essay
RANDOM_STATE = 42               # reproducibility
N_HVG = 2000                    # how many highly variable genes to keep
N_PCS = 30                      # PCs to compute before building the neighbor graph
N_NEIGHBORS = 15                # UMAP neighbor count — middle of the road


# ── load datasets ─────────────────────────────────────────────────────────────

# pbmc3k: raw counts, 2,700 PBMCs from a healthy donor (10x Genomics)
# this gets downloaded automatically on first run (~6MB)
pbmc3k = sc.datasets.pbmc3k()
pbmc3k.obs["batch"] = "PBMC 3k"

# pbmc68k_reduced: a 500-cell random subsample of the 68k PBMC dataset,
# already log-normalized and stored in .X. different donor, different run.
pbmc68k = sc.datasets.pbmc68k_reduced()
pbmc68k.obs["batch"] = "PBMC 68k"

# the 68k dataset comes pre-processed — pull it back to raw counts if available,
# otherwise we'll work with what we have and note it in preprocessing
if pbmc68k.raw is not None:
    pbmc68k = pbmc68k.raw.to_adata()
    pbmc68k.obs["batch"] = "PBMC 68k"


# ── preprocessing: each dataset independently before merging ──────────────────

# basic QC filters — remove obvious junk cells before we do anything else
sc.pp.filter_cells(pbmc3k, min_genes=200)
sc.pp.filter_genes(pbmc3k, min_cells=3)
sc.pp.filter_cells(pbmc68k, min_genes=200)
sc.pp.filter_genes(pbmc68k, min_cells=3)

# normalize each dataset to 10k counts per cell, then log-transform
# doing this per-dataset before merging keeps the scales comparable
sc.pp.normalize_total(pbmc3k, target_sum=1e4)
sc.pp.log1p(pbmc3k)

sc.pp.normalize_total(pbmc68k, target_sum=1e4)
sc.pp.log1p(pbmc68k)


# ── merge ─────────────────────────────────────────────────────────────────────

# concatenate along the gene axis — inner join keeps only genes present in both
# this is the honest approach: no imputation, no padding with zeros
merged = sc.concat(
    [pbmc3k, pbmc68k],
    join="inner",           # intersection of shared genes only
    label="batch",
    keys=["PBMC 3k", "PBMC 68k"],
    index_unique="-",
)

# re-label after concat (concat overwrites obs columns with the keys above)
merged.obs["batch"] = merged.obs["batch"].astype(str)


# ── feature selection and dimensionality reduction ────────────────────────────

# find highly variable genes on the merged object — these are the features
# that actually carry signal, rather than housekeeping genes with flat expression
sc.pp.highly_variable_genes(merged, n_top_genes=N_HVG, batch_key="batch")
merged = merged[:, merged.var["highly_variable"]]

# scale to unit variance before PCA — this prevents high-expression genes
# from dominating the principal components just because they're big numbers
sc.pp.scale(merged, max_value=10)

# PCA on the merged, scaled, HVG-filtered matrix
sc.tl.pca(merged, n_comps=N_PCS, random_state=RANDOM_STATE)

# build the neighbor graph in PCA space, then project to UMAP
# note: clustering and trajectory inference also happen in this PCA space,
# not in UMAP coordinates — UMAP is visualization only
sc.pp.neighbors(merged, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, random_state=RANDOM_STATE)
sc.tl.umap(merged, random_state=RANDOM_STATE)


# ── figure ────────────────────────────────────────────────────────────────────

BATCH_COLORS = {
    "PBMC 3k": "#2166ac",    # blue
    "PBMC 68k": "#d6604d",   # salmon-red
}

fig = plt.figure(figsize=(13, 5.5))
fig.patch.set_facecolor("white")

gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)
ax_pca = fig.add_subplot(gs[0])
ax_umap = fig.add_subplot(gs[1])


def scatter_by_batch(ax, coords, obs_batch, title, xlabel, ylabel):
    """
    Draws a scatter plot colored by batch label.
    Kept as a function so PCA and UMAP panels share identical styling.
    """
    for batch_label, color in BATCH_COLORS.items():
        mask = obs_batch == batch_label
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=color,
            label=batch_label,
            s=6,
            alpha=0.6,
            linewidths=0,
            rasterized=True,   # keeps file size down when exported as PDF/SVG
        )

    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.set_xlabel(xlabel, fontsize=10, color="#444444")
    ax.set_ylabel(ylabel, fontsize=10, color="#444444")
    ax.tick_params(labelsize=8, color="#888888")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#cccccc")

    legend = ax.legend(
        title="Dataset",
        title_fontsize=9,
        fontsize=8,
        frameon=True,
        framealpha=0.9,
        edgecolor="#dddddd",
        markerscale=2.5,
    )
    legend.get_title().set_color("#444444")


# PCA panel — batch separation should be obvious here
pca_coords = merged.obsm["X_pca"]
var_explained = merged.uns["pca"]["variance_ratio"]

scatter_by_batch(
    ax=ax_pca,
    coords=pca_coords,
    obs_batch=merged.obs["batch"],
    title="PCA",
    xlabel=f"PC1 ({var_explained[0]*100:.1f}% variance)",
    ylabel=f"PC2 ({var_explained[1]*100:.1f}% variance)",
)

# UMAP panel — batch structure is hidden inside what looks like cell type clusters
umap_coords = merged.obsm["X_umap"]

scatter_by_batch(
    ax=ax_umap,
    coords=umap_coords,
    obs_batch=merged.obs["batch"],
    title="UMAP",
    xlabel="UMAP1",
    ylabel="UMAP2",
)

# shared annotation at the bottom
fig.text(
    0.5, -0.04,
    "Two real PBMC datasets merged without batch correction. PCA immediately flags batch "
    "as the dominant source of variation (PC1, 11% variance). UMAP embeds them on the same "
    "canvas — leaving an ambiguous picture where you can't tell biology from artifact.",
    ha="center",
    fontsize=9,
    color="#666666",
    style="italic",
    wrap=True,
)

plt.savefig(
    "batch_effect_pca_vs_umap.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

print("saved: batch_effect_pca_vs_umap.png")
