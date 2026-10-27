# Scripts

| Script | Purpose | Runs on a fresh clone? |
|---|---|---|
| `reproduce_volcano.py` | Re-creates the volcano plot from `data/retinoblastoma_gene_results.csv` | **Yes** (tested) |
| `enrichment_analysis.py` | Gene-symbol annotation, over-representation analysis and preranked GSEA | No. Needs MSigDB gene-set tables and edited paths |
| `make_figures.py` | Volcano with symbols, MA plot, enrichment summary, interactive Plotly volcano | No. Paths point to the original environment |
| `make_extra_figures.py` | Enrichment heatmap, dot plots, GSEA curves, pathway gene heatmap, QC panels | No. Same reason |
| `build_database.py` | Builds `data/retinoblastoma.db` from the results tables | No. Same reason |
| `pca_heatmap_template.py` | PCA and sample heatmap for a raw count matrix | Template, tested on dummy data only |

The non-reproducible scripts were run in a temporary environment where the MSigDB v7.5.1 gene-set tables (Hallmark, KEGG, Reactome, GO BP) were loaded from local files. To redo the enrichment yourself, get the same gene sets from the R package `msigdbr` or from `gseapy`, then edit the file paths at the top of each script.
