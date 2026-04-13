# Dimensionality Reduction, Part 4: Which Method Do You Actually Use?

Code to reproduce the batch effect figure from [Part 4 of the Practical Guide to Dimensionality Reduction](https://quantavo.ai/blog/dimensionality-reduction-part-4) on the Quantavo blog.

## What this produces

A side-by-side PCA vs. UMAP plot using two real PBMC datasets (3k and 68k, downsampled) merged without batch correction. The figure demonstrates how PCA immediately surfaces batch as the dominant source of variation while UMAP buries it inside what appears to be clean cell-type clusters.

## Usage

```bash
pip install -r requirements.txt
python part4_batch_effect_demo.py
```

Saves `batch_effect_pca_vs_umap.png` to the current directory.

The datasets (pbmc3k and pbmc68k_reduced) are downloaded automatically by scanpy on first run.

## Dependencies

- `scanpy` — single-cell analysis and built-in PBMC datasets
- `matplotlib` — figure generation
- `numpy` — array operations

## Notes

- `RANDOM_STATE = 42` for reproducibility
- `N_HVG = 2000` highly variable genes selected after merging
- `N_PCS = 30` principal components for neighbor graph construction
- `N_NEIGHBORS = 15` for UMAP

## Part of the series

- [Part 1](https://quantavo.ai/blog/dimensionality-reduction-part-1) — Why your data has too many dimensions
- [Part 2](https://quantavo.ai/blog/dimensionality-reduction-part-2) — How PCA works
- [Part 3](https://quantavo.ai/blog/dimensionality-reduction-part-3) — t-SNE and UMAP
- [Part 4](https://quantavo.ai/blog/dimensionality-reduction-part-4) — Which method do you actually use?
