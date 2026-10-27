import pandas as pd, numpy as np, textwrap, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
import gseapy as gp
OUT='/mnt/user-data/outputs/'
d=pd.read_csv('/home/claude/annotated.csv')
ora=pd.read_csv(OUT+'enrichment_ORA.csv'); gsea=pd.read_csv(OUT+'enrichment_GSEA.csv')

def _fix(t):
    for a,b in [('Dna','DNA'),(' dna',' DNA'),('E2f','E2F'),('G2m','G2M'),('Myc','MYC'),('Tnfa','TNFα'),('Nfkb','NF-κB'),('nfkb','NF-κB'),('Kras','KRAS'),('Uv ','UV '),('Tca','TCA'),('Hdr','HDR'),('hrr','HRR'),('Slc ','SLC '),('Rna','RNA'),('Mtorc1','mTORC1'),('Hcm','HCM'),('Arvc','ARVC'),('Atp','ATP'),('p53','p53'),('P53','p53')]:
        t=t.replace(a,b)
    return t
def nice(s,w=40):
    for p in ['HALLMARK_','KEGG_','REACTOME_','GOBP_']: s=s.replace(p,'')
    return _fix(textwrap.fill(s.replace('_',' ').lower().capitalize(),w))
RED,BLUE='#d1495b','#2e6fbb'

# ---------- 1. Enrichment heatmap (Hallmark + KEGG, terms x direction)
def sgn(r): return (1 if r.direction=='Up' else -1)*-np.log10(max(r.fdr,1e-300))
rows=[]
for coll in ['Hallmark','KEGG','Reactome']:
    t=ora[(ora.collection==coll)&(ora.fdr<.05)]
    keep=pd.concat([t[t.direction=='Up'].nsmallest(10 if coll!='Reactome' else 8,'pvalue'),t[t.direction=='Down'].nsmallest(10 if coll!='Reactome' else 8,'pvalue')]).term.unique()
    for term in keep:
        r={'term':term,'coll':coll}
        for dirn in ['Up','Down']:
            x=ora[(ora.collection==coll)&(ora.term==term)&(ora.direction==dirn)]
            r[dirn]=(-np.log10(x.fdr.iloc[0]) if len(x) else 0)*(1 if dirn=='Up' else -1)
        rows.append(r)
H=pd.DataFrame(rows); H['net']=H.Up+H.Down
H=H.sort_values(['coll','net'],ascending=[True,False]); 
fig,axs=plt.subplots(1,3,figsize=(15,9),dpi=200)
for ax,coll in zip(axs,['Hallmark','KEGG','Reactome']):
    h=H[H.coll==coll].set_index('term')[['Up','Down']]; h.index=[nice(i,38) for i in h.index]
    sns.heatmap(h,cmap='vlag',center=0,vmin=-12,vmax=12,ax=ax,cbar=ax is axs[2],linewidths=.5,linecolor='white',
                cbar_kws={'label':'signed −log$_{10}$ FDR (capped at 12)\nred = enriched in up list, blue = in down list','shrink':.5} if ax is axs[2] else None)
    ax.set_title(coll+' (ORA)',loc='left'); ax.set_xticklabels(['Upregulated\ngenes','Downregulated\ngenes'],rotation=0); ax.tick_params(axis='y',labelsize=7.5)
plt.tight_layout(); plt.savefig(OUT+'enrichment_heatmap.png'); plt.close()

# ---------- 2. Dot plots (GO BP & Reactome)
fig,axs=plt.subplots(2,2,figsize=(13,10),dpi=200)
for i,coll in enumerate(['GO_BP','Reactome']):
    for j,dirn in enumerate(['Up','Down']):
        t=ora[(ora.collection==coll)&(ora.direction==dirn)&(ora.fdr<.05)].nsmallest(12,'pvalue').iloc[::-1]
        ax=axs[i,j]; sc=ax.scatter(t.overlap/t.list_size*100,[nice(x,45) for x in t.term],s=t.overlap*3.2,c=-np.log10(t.fdr),cmap='Reds' if dirn=='Up' else 'Blues',edgecolor='k',linewidth=.4)
        ax.set_xlabel('% of input genes in term'); ax.set_title(f'{coll.replace("_"," ")}: {dirn}-regulated',loc='left',fontsize=11)
        ax.tick_params(axis='y',labelsize=7.5); ax.margins(x=.15,y=.06); plt.colorbar(sc,ax=ax,label='−log$_{10}$ FDR',shrink=.7)
        for sp in ['top','right']: ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'enrichment_dotplots.png'); plt.close()

# ---------- 3. GSEA running-score curves
gg=pd.read_pickle('/tmp/msig_msigdbr_geneset_genes.pkl'); gs=pd.read_pickle('/tmp/msig_msigdbr_genesets.pkl'); g=pd.read_pickle('/tmp/msig_msigdbr_genes.pkl')
sym=g.drop_duplicates('gene_id').set_index('gene_id')['human_gene_symbol']; gg['sym']=gg.gene_id.map(sym)
want={'HALLMARK_E2F_TARGETS':('H',''),'HALLMARK_G2M_CHECKPOINT':('H',''),'HALLMARK_MYC_TARGETS_V1':('H',''),
      'REACTOME_THE_PHOTOTRANSDUCTION_CASCADE':('C2','CP:REACTOME'),'HALLMARK_OXIDATIVE_PHOSPHORYLATION':('H',''),'KEGG_CELL_CYCLE':('C2','CP:KEGG')}
sets={}
for name,(c,sc) in want.items():
    gid=gs[(gs.gs_name==name)].gs_id.iloc[0]; sets[name]=sorted(set(gg[gg.gs_id==gid].sym.dropna()))
rk=d.dropna(subset=['symbol']).copy()
rk=rk.reindex(rk.stat.abs().sort_values(ascending=False).index).drop_duplicates('symbol').set_index('symbol').stat.sort_values(ascending=False)
pr=gp.prerank(rnk=rk,gene_sets=sets,min_size=15,max_size=1000,permutation_num=1000,seed=42,threads=4,outdir=None,verbose=False)
fig,axs=plt.subplots(2,3,figsize=(14,7),dpi=200)
for ax,name in zip(axs.flat,sets):
    r=pr.results[name]; res=np.array(r['RES']); hits=np.array(r['hits']); n=len(res)
    pos=r['nes']>0; c=RED if pos else BLUE
    ax.plot(np.arange(n),res,color=c,lw=1.8); ax.axhline(0,c='k',lw=.6)
    ymin=min(res.min(),0); ax.vlines(hits,ymin-.07*(res.max()-ymin)-0.0,ymin-.01,color='k',lw=.4)
    ax.set_title(f'{nice(name,34)}\nNES {r["nes"]:.2f}, FDR {"< 0.001" if r["fdr"]<0.001 else format(r["fdr"],".3f")}',fontsize=9.5,loc='left')
    ax.set_xlabel('Gene rank (high → low statistic)',fontsize=8); ax.set_ylabel('Running enrichment score',fontsize=8)
    for sp in ['top','right']: ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'gsea_running_scores.png'); plt.close()

# ---------- 4. Gene x pathway heatmap coloured by log2FC
pw=['HALLMARK_E2F_TARGETS','HALLMARK_G2M_CHECKPOINT','HALLMARK_MYC_TARGETS_V1','KEGG_CELL_CYCLE','REACTOME_THE_PHOTOTRANSDUCTION_CASCADE','HALLMARK_OXIDATIVE_PHOSPHORYLATION']
lfc=d.dropna(subset=['symbol']).drop_duplicates('symbol').set_index('symbol')
sigl=lfc[(lfc.padj<.05)&(lfc.log2foldchange.abs()>=1)]
genes={}
for name in pw[:4]+pw[5:]:
    s=sigl[sigl.index.isin(sets[name])].log2foldchange
    genes[name]=s.sort_values(ascending=(name==pw[5])).index[:12 if name!=pw[3] else 8]
photo=sigl[sigl.index.isin(sets[pw[4]])].log2foldchange.sort_values().index[:12]; genes[pw[4]]=photo
sel=list(dict.fromkeys([x for k in genes.values() for x in k]))
M=pd.DataFrame(np.nan,index=sel,columns=pw)
for name in pw:
    for gname in sel:
        if gname in sets[name] and gname in lfc.index: M.loc[gname,name]=lfc.loc[gname,'log2foldchange']
M=M.loc[M.notna().any(axis=1)]
order=M.apply(lambda r: np.nanmean(r),axis=1).sort_values(ascending=False).index; M=M.loc[order]
fig,ax=plt.subplots(figsize=(7.5,11),dpi=200)
sns.heatmap(M,cmap='vlag',center=0,ax=ax,linewidths=.4,linecolor='white',cbar_kws={'label':'log$_2$ fold change','shrink':.5},mask=M.isna(),yticklabels=True)
ax.set_facecolor('#f1f1f1'); ax.set_xticklabels([nice(c,18) for c in pw],rotation=40,ha='right',fontsize=8); ax.tick_params(axis='y',labelsize=7)
ax.set_title('Top differentially expressed genes in key pathways\n(grey = gene not in pathway; all genes padj < 0.05)',fontsize=10,loc='left')
plt.tight_layout(); plt.savefig(OUT+'pathway_gene_heatmap.png'); plt.close()

# ---------- 5. Top genes diverging bars (by padj) + QC panels
sig=d[d.status.isin(['Up','Down'])]
top=pd.concat([sig[sig.status=='Up'].nsmallest(15,'padj'),sig[sig.status=='Down'].nsmallest(15,'padj')]).sort_values('log2foldchange')
fig,axs=plt.subplots(1,3,figsize=(16,6.5),dpi=200,gridspec_kw={'width_ratios':[1.2,1,1]})
ax=axs[0]; ax.barh(top.label,top.log2foldchange,color=[RED if v>0 else BLUE for v in top.log2foldchange])
ax.axvline(0,c='k',lw=.6); ax.set_xlabel('log$_2$ fold change'); ax.set_title('A  Top 15 up / down genes (by padj)',loc='left',fontsize=11); ax.tick_params(axis='y',labelsize=8)
raw=pd.read_csv('/mnt/user-data/uploads/1793077542523_retinoblastoma_gene_results.csv')
ax=axs[1]; ax.hist(raw.pvalue.dropna(),bins=50,color='#6b7280'); ax.set_xlabel('p-value'); ax.set_ylabel('Genes'); ax.set_title('B  P-value histogram (QC)',loc='left',fontsize=11)
ax=axs[2]; cnt=d.status.value_counts().reindex(['Up','Down','Sig. (small FC)','NS']); cc={'Up':RED,'Down':BLUE,'Sig. (small FC)':'#8fb8de','NS':'#c9ccd1'}
ax.bar(cnt.index,cnt.values,color=[cc[i] for i in cnt.index]); ax.set_yscale('log'); ax.set_ylabel('Genes (log scale)'); ax.set_title('C  Genes per category',loc='left',fontsize=11)
for i,v in enumerate(cnt.values): ax.text(i,v*1.1,f'{v:,}',ha='center',fontsize=9)
ax.tick_params(axis='x',labelsize=8)
for a in axs:
    for sp in ['top','right']: a.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'top_genes_and_qc.png'); plt.close()
print('ok'); print(M.shape, {k:len(v) for k,v in genes.items()})
