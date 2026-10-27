import pandas as pd, numpy as np, textwrap
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text
import plotly.graph_objects as go
OUT='/mnt/user-data/outputs/'
d=pd.read_csv('/home/claude/annotated.csv')
d['nl']=-np.log10(d.padj.clip(lower=1e-300))
col={'NS':'#c9ccd1','Sig. (small FC)':'#8fb8de','Up':'#d1495b','Down':'#2e6fbb'}
known=['RB1','MYCN','NRL','GRK1','DLL3','MKI67']
print(d[d.label.isin(known)][['label','log2foldchange','padj','status']])

# ---- volcano with symbols
fig,ax=plt.subplots(figsize=(10,8),dpi=200)
for s in ['NS','Sig. (small FC)','Down','Up']:
    x=d[d.status==s]; ax.scatter(x.log2foldchange,x.nl,s=9 if s!='NS' else 6,c=col[s],alpha=.75 if s!='NS' else .5,linewidths=0,label=f'{s} (n={len(x):,})')
for v in [1,-1]: ax.axvline(v,ls='--',c='#555',lw=.9)
ax.axhline(-np.log10(.05),ls='--',c='#555',lw=.9)
sig=d[d.status.isin(['Up','Down'])]
lab=pd.concat([sig[sig.status=='Up'].nsmallest(10,'padj'),sig[sig.status=='Down'].nsmallest(10,'padj'),d[d.label.isin(['RB1','MYCN'])]]).drop_duplicates('gene')
t=[ax.text(r.log2foldchange,r.nl,r.label,fontsize=8,fontweight='bold' if r.label in ('RB1','MYCN') else 'normal') for _,r in lab.iterrows()]
adjust_text(t,ax=ax,arrowprops=dict(arrowstyle='-',color='#777',lw=.5),expand=(1.3,1.6))
ax.set_xlabel('log$_2$ fold change',fontsize=12); ax.set_ylabel('−log$_{10}$ adjusted p-value',fontsize=12)
ax.set_title('Retinoblastoma differential expression (DESeq2)\npadj < 0.05 and |log$_2$FC| ≥ 1',fontsize=13)
ax.legend(frameon=False,loc='upper left',fontsize=9,markerscale=2)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'volcano_plot.png'); plt.savefig(OUT+'volcano_plot.pdf'); plt.close()

# ---- MA plot
fig,ax=plt.subplots(figsize=(8,5.5),dpi=200)
for s in ['NS','Sig. (small FC)','Down','Up']:
    x=d[d.status==s]; ax.scatter(np.log10(x.basemean+1),x.log2foldchange,s=5,c=col[s],alpha=.6,linewidths=0,label=s)
ax.axhline(0,c='k',lw=.6); ax.set_xlabel('log$_{10}$(mean normalised count + 1)'); ax.set_ylabel('log$_2$ fold change')
ax.set_title('MA plot'); ax.legend(frameon=False,markerscale=3,fontsize=8)
for sp in ['top','right']: ax.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'MA_plot.png'); plt.close()

# ---- enrichment figure
ora=pd.read_csv(OUT+'enrichment_ORA.csv'); gs=pd.read_csv(OUT+'enrichment_GSEA.csv')

def _fix(t):
    for a,b in [('Dna','DNA'),(' dna',' DNA'),('E2f','E2F'),('G2m','G2M'),('Myc','MYC'),('Tnfa','TNFα'),('Nfkb','NF-κB'),('nfkb','NF-κB'),('Kras','KRAS'),('Uv ','UV '),('Tca','TCA'),('Hdr','HDR'),('hrr','HRR'),('Slc ','SLC '),('Rna','RNA'),('Mtorc1','mTORC1'),('Hcm','HCM'),('Arvc','ARVC'),('Atp','ATP'),('p53','p53'),('P53','p53')]:
        t=t.replace(a,b)
    return t
def nice(s,pre): return _fix(textwrap.fill(s.replace(pre,'').replace('_',' ').lower().capitalize(),38))
def ora_panel(ax,coll,pre,n,title):
    rows=[]
    for dirn,sgn in [('Up',1),('Down',-1)]:
        t=ora[(ora.collection==coll)&(ora.direction==dirn)&(ora.fdr<.05)].nsmallest(n,'pvalue').copy()
        t['v']=sgn*-np.log10(t.fdr.clip(lower=1e-300)); rows.append(t)
    t=pd.concat(rows).sort_values('v')
    ax.barh([nice(x,pre) for x in t.term],t.v,color=['#2e6fbb' if v<0 else '#d1495b' for v in t.v])
    ax.axvline(0,c='k',lw=.7); ax.set_xlabel('Signed −log$_{10}$ FDR  (← down | up →)'); ax.set_title(title,fontsize=11,loc='left')
    ax.tick_params(axis='y',labelsize=7.5)
    for sp in ['top','right']: ax.spines[sp].set_visible(False)
fig=plt.figure(figsize=(15,9),dpi=200)
ax1=fig.add_subplot(1,3,1); ora_panel(ax1,'Hallmark','HALLMARK_',8,'A  Hallmark (ORA)')
ax2=fig.add_subplot(1,3,2); ora_panel(ax2,'GO_BP','GOBP_',8,'B  GO Biological Process (ORA)')
ax3=fig.add_subplot(1,3,3)
h=gs[(gs.collection=='Hallmark')&(gs['FDR q-val']<.05)].sort_values('NES')
h=pd.concat([h.head(8),h.tail(8)])
ax3.barh([nice(x,'HALLMARK_') for x in h.Term],h.NES,color=['#2e6fbb' if v<0 else '#d1495b' for v in h.NES])
ax3.axvline(0,c='k',lw=.7); ax3.set_xlabel('Normalised enrichment score (NES)'); ax3.set_title('C  Hallmark (preranked GSEA)',fontsize=11,loc='left'); ax3.tick_params(axis='y',labelsize=7.5)
for sp in ['top','right']: ax3.spines[sp].set_visible(False)
plt.tight_layout(); plt.savefig(OUT+'enrichment_plots.png'); plt.close()

# ---- interactive plotly volcano
fig=go.Figure()
order=['NS','Sig. (small FC)','Down','Up']
for s in order:
    x=d[d.status==s]
    fig.add_trace(go.Scattergl(x=x.log2foldchange,y=x.nl,mode='markers',name=f'{s} ({len(x):,})',
        marker=dict(size=5 if s!='NS' else 4,color=col[s],opacity=.8 if s!='NS' else .45),
        customdata=np.c_[x.label,x.gene,x.basemean.round(1),x.padj],
        hovertemplate='<b>%{customdata[0]}</b><br>%{customdata[1]}<br>log2FC: %{x:.2f}<br>padj: %{customdata[3]:.2e}<br>baseMean: %{customdata[2]}<extra></extra>'))
ann=[dict(x=r.log2foldchange,y=r.nl,text=r.label,showarrow=True,arrowhead=0,arrowwidth=.7,ax=18,ay=-18,font=dict(size=10)) for _,r in lab.head(0).iterrows()]
top=pd.concat([sig[sig.status=='Up'].nsmallest(6,'padj'),sig[sig.status=='Down'].nsmallest(6,'padj'),d[d.label.isin(['RB1','MYCN'])]]).drop_duplicates('gene')
ann=[dict(x=r.log2foldchange,y=r.nl,text=r.label,showarrow=True,arrowhead=0,arrowwidth=.7,ax=22,ay=-16,font=dict(size=10)) for _,r in top.iterrows()]
fig.update_layout(template='plotly_white',title='Retinoblastoma differential expression — interactive volcano plot',
    xaxis_title='log₂ fold change',yaxis_title='−log₁₀ adjusted p-value',height=700,annotations=ann,
    shapes=[dict(type='line',x0=v,x1=v,yref='paper',y0=0,y1=1,line=dict(dash='dash',color='#666',width=1)) for v in (-1,1)]+
           [dict(type='line',xref='paper',x0=0,x1=1,y0=-np.log10(.05),y1=-np.log10(.05),line=dict(dash='dash',color='#666',width=1))],
    legend=dict(title='Click to show/hide',x=.01,y=.99),margin=dict(t=70))
fig.write_html(OUT+'volcano_interactive.html',include_plotlyjs=True,config=dict(displaylogo=False,scrollZoom=True))
print('done')
