"""PCA + heatmap from a raw count matrix.  Usage:
   python pca_heatmap_template.py counts.csv metadata.csv significant_DEGs_annotated.csv
 counts.csv   : rows = Ensembl gene IDs (first column), columns = samples (raw counts)
 metadata.csv : columns 'sample','condition'  (e.g. tumour / normal)
"""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, seaborn as sns
from sklearn.decomposition import PCA
counts=pd.read_csv(sys.argv[1],index_col=0); meta=pd.read_csv(sys.argv[2]).set_index('sample').loc[counts.columns]
deg=pd.read_csv(sys.argv[3])
counts=counts[counts.sum(axis=1)>10]
logcpm=np.log2(counts/counts.sum()*1e6+1)            # simple stand-in for DESeq2 VST
# PCA on the 2,000 most variable genes
top=logcpm.var(axis=1).nlargest(2000).index
p=PCA(2).fit(logcpm.loc[top].T); pc=p.transform(logcpm.loc[top].T)
fig,ax=plt.subplots(figsize=(6,5),dpi=200)
for c,g in meta.groupby('condition'):
    i=[list(counts.columns).index(s) for s in g.index]; ax.scatter(pc[i,0],pc[i,1],s=60,label=c)
ax.set_xlabel(f'PC1 ({p.explained_variance_ratio_[0]*100:.1f}%)'); ax.set_ylabel(f'PC2 ({p.explained_variance_ratio_[1]*100:.1f}%)'); ax.legend(frameon=False); ax.set_title('PCA')
plt.tight_layout(); plt.savefig('pca_plot.png')
# heatmap: top 25 up + top 25 down by padj, row z-scores
sig=deg.sort_values('padj'); g=pd.concat([sig[sig.status=='Up'].head(25),sig[sig.status=='Down'].head(25)])
g=g[g.gene.isin(logcpm.index)]; m=logcpm.loc[g.gene]; m.index=g.symbol.fillna(g.gene).values
z=m.sub(m.mean(axis=1),axis=0).div(m.std(axis=1),axis=0)
lut=dict(zip(meta.condition.unique(),sns.color_palette('Set2')))
cg=sns.clustermap(z,cmap='vlag',center=0,col_colors=meta.condition.map(lut).values,figsize=(8,10),yticklabels=True,xticklabels=False)
cg.savefig('heatmap_top_genes.png',dpi=200)
