# Dimensionality Reduction, Part 1 — Companion Code

**Blog post:** [The Practical Guide to Dimensionality Reduction, Part 1: Why Your Data Has Too Many Dimensions](https://quantavo.ai/blog/dimensionality-reduction-part-1)

This script generates a two-panel figure comparing a raw two-gene projection of 600 simulated cells against the same data after PCA on highly variable genes, illustrating why dimensionality reduction changes what you can see in your data.

## How to run

```bash
pip install -r requirements.txt
python part1_create_dim_reduction_comparison.py
```

Output: `dim-reduction-comparison.png` in the current directory.

> **Note:** The figure may differ slightly from the one in the blog post due to differences in library versions or rendering environments. The underlying data and logic are identical.
