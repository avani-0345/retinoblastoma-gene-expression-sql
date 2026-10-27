"""Re-create the volcano plot from the input DESeq2 table.

Run from the repository root:
    pip install -r requirements.txt
    python scripts/reproduce_volcano.py

Input : data/retinoblastoma_gene_results.csv   (DESeq2 results, Ensembl IDs)
Output: figures/volcano_plot_reproduced.png
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PADJ, LFC = 0.05, 1.0

raw = pd.read_csv(ROOT / "data" / "retinoblastoma_gene_results.csv")
df = raw.dropna(subset=["padj", "log2foldchange"]).copy()
df["nl"] = -np.log10(df.padj.clip(lower=1e-300))

def category(r):
    if r.padj < PADJ and r.log2foldchange >= LFC: return "Up"
    if r.padj < PADJ and r.log2foldchange <= -LFC: return "Down"
    if r.padj < PADJ: return "Sig. (small FC)"
    return "NS"

df["status"] = df.apply(category, axis=1)
colors = {"NS": "#c9ccd1", "Sig. (small FC)": "#8fb8de", "Up": "#d1495b", "Down": "#2e6fbb"}

fig, ax = plt.subplots(figsize=(10, 8), dpi=200)
for s in ["NS", "Sig. (small FC)", "Down", "Up"]:
    x = df[df.status == s]
    ax.scatter(x.log2foldchange, x.nl, s=8, c=colors[s], alpha=0.7, linewidths=0, label=f"{s} (n={len(x):,})")
for v in (-LFC, LFC): ax.axvline(v, ls="--", c="#555", lw=0.9)
ax.axhline(-np.log10(PADJ), ls="--", c="#555", lw=0.9)
ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("Retinoblastoma differential expression (DESeq2)")
ax.legend(frameon=False, loc="upper left", markerscale=2)
plt.tight_layout()
out = ROOT / "figures" / "volcano_plot_reproduced.png"
plt.savefig(out); print("Saved", out); print(df.status.value_counts().to_string())
