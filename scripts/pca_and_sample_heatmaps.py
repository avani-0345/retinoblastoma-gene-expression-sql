"""PCA and sample-level heatmaps from the GEO raw counts (GSE125903).

Run from the repository root:   python scripts/pca_and_sample_heatmaps.py
Inputs : data/GSE125903_raw_counts_GRCh38_p13_NCBI.tsv.gz  (NCBI-generated raw counts, Entrez IDs, 7 samples)
         data/entrez_to_symbol.csv                         (Entrez ID -> gene symbol)
         data/annotated_results_all_genes.csv              (DESeq2 table with symbols and categories)
Outputs: figures/pca_samples.png, figures/heatmap_key_genes.png, figures/heatmap_top_DE_genes_samples.png

NOTE: sample labels (tumour vs retina) are NOT included in the counts file. Two samples (GSM3584494, GSM3584495)
show a photoreceptor-high, low-proliferation profile and separate from the other five on PC1. Samples are therefore
grouped by their expression pattern (k-means on PC1/PC2), not by official GEO labels. Check the GEO sample list.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
cnt = pd.read_csv(ROOT/"data"/"GSE125903_raw_counts_GRCh38_p13_NCBI.tsv.gz", sep="\t", index_col=0)
sym = pd.read_csv(ROOT/"data"/"entrez_to_symbol.csv").drop_duplicates("GeneID").set_index("GeneID").symbol
de  = pd.read_csv(ROOT/"data"/"annotated_results_all_genes.csv")

lib = cnt.sum()
cpm = cnt / lib * 1e6                                  # counts per million, all genes in library size
logcpm = np.log2(cpm + 1)
logcpm["symbol"] = sym.reindex(logcpm.index)
logcpm = logcpm.dropna(subset=["symbol"]).groupby("symbol").mean()   # one row per symbol (duplicates averaged)
samples = list(cnt.columns)
short = {s: "GSM..." + s[-4:] for s in samples}
lab = [short[s] for s in samples]

# ---- PCA: 2,000 most variable expressed genes
expr = logcpm[(cpm.groupby(sym.reindex(cpm.index)).mean() > 1).reindex(logcpm.index).any(axis=1)] if False else logcpm
keep = logcpm[(logcpm > np.log2(2)).sum(axis=1) >= 3]          # CPM >= 1 in at least 3 samples
top = keep.var(axis=1).nlargest(2000).index
p = PCA(3).fit(keep.loc[top].T.values); pc = p.transform(keep.loc[top].T.values)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.8), dpi=200)
km = KMeans(2, n_init=10, random_state=0).fit(pc[:, :2]); grp = km.labels_
mk = logcpm.loc["MKI67", samples].values
hi = int(np.argmax([mk[grp == k].mean() for k in (0, 1)]))
gname = {hi: "Proliferative profile (high MKI67)", 1 - hi: "Photoreceptor-high, low proliferation"}
gcol = {hi: "#d1495b", 1 - hi: "#2e6fbb"}
for k in (0, 1):
    m = grp == k; ax[0].scatter(pc[m, 0], pc[m, 1], s=120, c=gcol[k], edgecolor="k", zorder=3, label=f"{gname[k]} (n={m.sum()})")
for i, s_ in enumerate(samples): ax[0].annotate(short[s_], (pc[i, 0], pc[i, 1]), xytext=(7, 6), textcoords="offset points", fontsize=8)
ax[0].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.14))
ax[0].set_xlabel(f"PC1 ({p.explained_variance_ratio_[0]*100:.1f}%)"); ax[0].set_ylabel(f"PC2 ({p.explained_variance_ratio_[1]*100:.1f}%)")
ax[0].set_title("PCA of 7 samples\n(top 2,000 variable genes, log2 CPM)", fontsize=11, loc="left"); ax[0].axhline(0, c="#ccc", lw=.6); ax[0].axvline(0, c="#ccc", lw=.6)
ax[1].bar(lab, lib / 1e6, color=[gcol[g] for g in grp]); ax[1].set_ylabel("Library size (million reads counted)"); ax[1].tick_params(axis="x", labelsize=7, rotation=30); ax[1].set_title("Sequencing depth per sample", fontsize=11, loc="left")
for a in ax:
    for sp in ["top", "right"]: a.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(ROOT/"figures"/"pca_samples.png"); plt.close()

# ---- Heatmap 1: key genes (absolute log2 CPM)
groups = {"Proliferation / E2F targets": ["MKI67", "TOP2A", "E2F1", "MCM2", "CDC20", "NDC80"],
          "RB1 / MYCN / Notch": ["RB1", "MYCN", "DLL3"],
          "Photoreceptor / phototransduction": ["CRX", "NRL", "RHO", "SAG", "GNAT1", "PDE6B", "GRK1", "RCVRN"]}
genes = [g for v in groups.values() for g in v if g in logcpm.index]
M = logcpm.loc[genes, samples]; M.columns = lab
fig, ax = plt.subplots(figsize=(7.5, 8), dpi=200)
sns.heatmap(M, cmap="YlOrRd", annot=True, fmt=".1f", annot_kws={"size": 8}, linewidths=.5, linecolor="white", ax=ax, cbar_kws={"label": "log2 (CPM + 1)", "shrink": .6})
y = 0
for name, gl in groups.items():
    n = len([g for g in gl if g in logcpm.index]); ax.hlines(y, -.8, 0, color="k", lw=0)
    ax.text(len(samples) + 2.3, y + n / 2, name, rotation=270, va="center", ha="center", fontsize=8.5, color="#14213d", fontweight="bold")
    y += n
    if y < len(genes): ax.axhline(y, color="k", lw=1.5)
ax.set_title("Expression of key genes across the 7 samples", fontsize=11, loc="left"); ax.tick_params(axis="y", labelsize=9, rotation=0); ax.tick_params(axis="x", labelsize=8)
ax.tick_params(axis="x", labelsize=7, rotation=30)
plt.tight_layout(); plt.savefig(ROOT/"figures"/"heatmap_key_genes.png"); plt.close()

# ---- Heatmap 2: top DE genes (row z-score across tumours)
sig = de[de.status.isin(["Up", "Down"])].dropna(subset=["symbol"]).drop_duplicates("symbol")
sel = pd.concat([sig[sig.status == "Up"].nsmallest(20, "padj"), sig[sig.status == "Down"].nsmallest(20, "padj")])
sel = sel[sel.symbol.isin(logcpm.index)]
Z = logcpm.loc[sel.symbol, samples]; Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1).replace(0, np.nan), axis=0); Z.columns = lab
rowc = sel.set_index("symbol").status.map({"Up": "#d1495b", "Down": "#2e6fbb"})
colc = pd.Series([gcol[g] for g in grp], index=Z.columns)
cg = sns.clustermap(Z, cmap="vlag", center=0, row_cluster=False, col_cluster=True, row_colors=rowc.loc[Z.index].values, col_colors=colc.values, figsize=(7.5, 11), yticklabels=True, cbar_kws={"label": "row z-score"})
cg.ax_heatmap.tick_params(axis="y", labelsize=8); cg.ax_heatmap.tick_params(axis="x", labelsize=7)
cg.fig.suptitle("Top 20 up (red bar) and top 20 down (blue bar) DE genes\nZ-scored across the 7 samples. Column bar: red = proliferative profile, blue = photoreceptor-high group", fontsize=10, y=1.03, x=0.58)
cg.savefig(ROOT/"figures"/"heatmap_top_DE_genes_samples.png", dpi=200, bbox_inches="tight"); plt.close("all")
print({short[x]: gname[g] for x, g in zip(samples, grp)}); print("PC1/PC2/PC3 variance:", np.round(p.explained_variance_ratio_ * 100, 1)); print("library sizes (M):", (lib / 1e6).round(1).to_dict())
print("genes plotted:", len(genes), "| DE genes in heatmap 2:", len(Z))
print(M.round(1).to_string())
