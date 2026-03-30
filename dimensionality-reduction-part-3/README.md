# Dimensionality Reduction Part 3 — t-SNE vs UMAP

Companion code for **[The Practical Guide to Dimensionality Reduction, Part 3: t-SNE and UMAP: When Linear Methods Aren't Enough](https://quantavo.ai/blog/dimensionality-reduction-part-3)** by Cory Henn.

## What this script does

Generates a side-by-side comparison of t-SNE and UMAP embeddings using the 10x Genomics PBMC 3k dataset. Both embeddings are computed from the same top 30 of 50 PCs for a fair, apples-to-apples comparison.

## Requirements

```
pip install -r requirements.txt
```

## Usage

```bash
python part3_tsne_umap_comparison.py
```

Outputs `tsne-vs-umap-comparison.png` in the current directory.

## Dataset

The 10x Genomics PBMC 3k dataset is downloaded automatically by `scanpy` on first run (~30 MB). No manual download required.

## Notes

- `random_state=42` is set throughout for reproducibility
- Both embeddings use `n_pcs=30` from a 50-component PCA
- t-SNE uses `perplexity=30` (scanpy default)
- Cell type labels are from the Louvain clustering in the preprocessed dataset
