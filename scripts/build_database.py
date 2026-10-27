import pandas as pd, sqlite3, os
OUT='/mnt/user-data/outputs/'
db=OUT+'retinoblastoma.db'
if os.path.exists(db): os.remove(db)
con=sqlite3.connect(db)
d=pd.read_csv('/home/claude/annotated.csv').rename(columns={'gene':'ensembl_id','basemean':'base_mean','log2foldchange':'log2fc'})
d[['ensembl_id','symbol','base_mean','log2fc','lfcse','stat','pvalue','padj','status']].to_sql('genes',con,index=False)
ora=pd.read_csv(OUT+'enrichment_ORA.csv'); ora.to_sql('enrichment_ora',con,index=False)
g=pd.read_csv(OUT+'enrichment_GSEA.csv').rename(columns={'Term':'term','NES':'nes','FDR q-val':'fdr','NOM p-val':'pvalue','ES':'es'})
g[['collection','term','es','nes','pvalue','fdr']].to_sql('enrichment_gsea',con,index=False)
# Hallmark membership
gs=pd.read_pickle('/tmp/msig_msigdbr_genesets.pkl'); gg=pd.read_pickle('/tmp/msig_msigdbr_geneset_genes.pkl'); gn=pd.read_pickle('/tmp/msig_msigdbr_genes.pkl')
sym=gn.drop_duplicates('gene_id').set_index('gene_id').human_gene_symbol
h=gs[gs.gs_cat=='H'][['gs_id','gs_name']]
m=gg[gg.gs_id.isin(h.gs_id)].merge(h,on='gs_id'); m['symbol']=m.gene_id.map(sym)
m[['gs_name','symbol']].dropna().drop_duplicates().rename(columns={'gs_name':'pathway'}).to_sql('hallmark_membership',con,index=False)
con.executescript('CREATE INDEX ix_g_sym ON genes(symbol); CREATE INDEX ix_g_status ON genes(status); CREATE INDEX ix_m_sym ON hallmark_membership(symbol);')
con.commit(); print(pd.read_sql("select name from sqlite_master where type='table'",con))
