import pandas as pd, numpy as np, warnings
from scipy.stats import hypergeom, false_discovery_control
warnings.filterwarnings('ignore')
SRC='/mnt/user-data/uploads/1793077542523_retinoblastoma_gene_results.csv'
OUT='/mnt/user-data/outputs/'
raw=pd.read_csv(SRC)

# --- 1. Ensembl -> symbol
g=pd.read_pickle('/tmp/msig_msigdbr_genes.pkl')
m=g.dropna(subset=['human_ensembl_gene']).drop_duplicates('human_ensembl_gene').set_index('human_ensembl_gene')['human_gene_symbol']
raw['symbol']=raw.gene.map(m)
d=raw.dropna(subset=['padj','log2foldchange']).copy()
print('symbol coverage: all',raw.symbol.notna().mean().round(3),'tested',d.symbol.notna().mean().round(3))
d['label']=d.symbol.fillna(d.gene)
def cls(r):
    if r.padj<0.05 and r.log2foldchange>=1: return 'Up'
    if r.padj<0.05 and r.log2foldchange<=-1: return 'Down'
    if r.padj<0.05: return 'Sig. (small FC)'
    return 'NS'
d['status']=d.apply(cls,axis=1)
d.to_csv('/home/claude/annotated.csv',index=False)
sig=d[d.status.isin(['Up','Down'])].sort_values('padj')
sig[['gene','symbol','basemean','log2foldchange','lfcse','pvalue','padj','status']].to_csv(OUT+'significant_DEGs_annotated.csv',index=False)
print('Top up:',sig[sig.status=='Up'].head(15)[['label','log2foldchange','padj']].to_string())
print('Top down:',sig[sig.status=='Down'].head(15)[['label','log2foldchange','padj']].to_string())

# --- 2. gene sets
gs=pd.read_pickle('/tmp/msig_msigdbr_genesets.pkl'); gg=pd.read_pickle('/tmp/msig_msigdbr_geneset_genes.pkl')
sym=g.drop_duplicates('gene_id').set_index('gene_id')['human_gene_symbol']
gg['sym']=gg.gene_id.map(sym); gg=gg.dropna(subset=['sym'])
cols={'Hallmark':('H',''),'KEGG':('C2','CP:KEGG'),'Reactome':('C2','CP:REACTOME'),'GO_BP':('C5','GO:BP')}
sets={}
for k,(c,sc) in cols.items():
    ids=gs[(gs.gs_cat==c)&(gs.gs_subcat==sc)].set_index('gs_id').gs_name
    x=gg[gg.gs_id.isin(ids.index)].groupby('gs_id').sym.apply(lambda s:set(s))
    sets[k]={ids[i]:v for i,v in x.items()}
    print(k,len(sets[k]))

# --- 3. ORA
bg=set(d.symbol.dropna())
lists={'Up':set(sig[sig.status=='Up'].symbol.dropna()),'Down':set(sig[sig.status=='Down'].symbol.dropna())}
print({k:len(v) for k,v in lists.items()},'background',len(bg))
rows=[]
for coll,S in sets.items():
    res=[]
    for name,genes in S.items():
        G=genes&bg
        if not 10<=len(G)<=500: continue
        for dirn,L in lists.items():
            k=len(L&G)
            if k<3: continue
            p=hypergeom.sf(k-1,len(bg),len(G),len(L))
            res.append((coll,dirn,name,k,len(G),len(L),k/len(G),p,';'.join(sorted(L&G)[:12])))
    r=pd.DataFrame(res,columns=['collection','direction','term','overlap','set_size','list_size','frac_of_set','pvalue','genes'])
    r['fdr']=r.groupby('direction').pvalue.transform(lambda p: false_discovery_control(p,method='bh'))
    rows.append(r)
ora=pd.concat(rows).sort_values(['collection','direction','pvalue'])
ora.to_csv(OUT+'enrichment_ORA.csv',index=False)
for coll in sets:
    for dirn in ['Up','Down']:
        t=ora[(ora.collection==coll)&(ora.direction==dirn)&(ora.fdr<0.05)]
        print(f'\n== {coll} {dirn}: {len(t)} terms FDR<0.05');print(t.head(8)[['term','overlap','set_size','fdr']].to_string(index=False))

# --- 4. preranked GSEA
import gseapy as gp
rk=d.dropna(subset=['symbol']).copy()
rk=rk.reindex(rk.stat.abs().sort_values(ascending=False).index).drop_duplicates('symbol').set_index('symbol').stat.sort_values(ascending=False)
print('ranked genes',len(rk))
allg=[]
for coll in ['Hallmark','KEGG','Reactome']:
    pr=gp.prerank(rnk=rk,gene_sets={k:sorted(v) for k,v in sets[coll].items()},min_size=15,max_size=500,permutation_num=1000,seed=42,threads=4,outdir=None,verbose=False)
    r=pr.res2d.copy(); r['collection']=coll; allg.append(r)
gsea=pd.concat(allg)
for c in ['ES','NES','NOM p-val','FDR q-val','FWER p-val']: gsea[c]=pd.to_numeric(gsea[c],errors='coerce')
gsea.to_csv(OUT+'enrichment_GSEA.csv',index=False)
for coll in ['Hallmark','KEGG','Reactome']:
    t=gsea[(gsea.collection==coll)&(gsea['FDR q-val']<0.05)].sort_values('NES')
    print(f'\n== GSEA {coll}: {len(t)} sig; pos NES {(t.NES>0).sum()}, neg {(t.NES<0).sum()}')
    print(t.head(6)[['Term','NES','FDR q-val']].to_string(index=False)); print(t.tail(6)[['Term','NES','FDR q-val']].to_string(index=False))
