# Dimensionality Reduction, Part 2 — Companion Code

**Blog post:** [The Practical Guide to Dimensionality Reduction, Part 2: PCA](https://quantavo.ai/blog/dimensionality-reduction-part-2)

This script generates a three-panel figure using the 10x Genomics PBMC 3k dataset (downloaded automatically by scanpy): an elbow plot showing variance explained per PC, a UMAP at 15 PCs (the correct elbow), and a UMAP at 100 PCs (too many), demonstrating how over-including PCs lets noise corrupt your embedding.

## How to run

```bash
pip install -r requirements.txt
python part2_elbow_pc_comparison.py
```

Output: `elbow-pc-comparison.png` in the current directory.

> **Note:** The figure may differ slightly from the one in the blog post due to differences in library versions or rendering environments. The underlying data and logic are identical.
