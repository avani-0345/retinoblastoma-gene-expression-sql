# Retinoblastoma RNA-seq: Differential Expression and Pathway Analysis

![Volcano plot](figures/volcano_plot.png)

Exploration of a DESeq2 differential expression table from public dataset GEO **GSE125903**, which compares retinoblastoma tumours with adult retina. Retinoblastoma is a childhood eye cancer caused by loss of the *RB1* gene. The analysis asks whether the data show the expected biology: **proliferation programmes up, retinal (photoreceptor) identity down**. They do.

**Contents:** [Key results](#key-results) · [Figures](#figures) · [SQL](#sql) · [Repository structure](#repository-structure) · [Reproducing](#reproducing) · [Limitations](#limitations) · [How this was made](#how-this-was-made)

📄 Full write-up: [`REPORT.pdf`](REPORT.pdf) · 🖼️ LinkedIn slides: [`docs/retinoblastoma_linkedin_carousel.pdf`](docs/retinoblastoma_linkedin_carousel.pdf)

## Key results
| | |
|---|---|
| Genes in input table | 57,773 |
| Genes with a valid adjusted p-value | 19,952 |
| **Upregulated** (padj < 0.05, log2FC ≥ 1) | **2,168** |
| **Downregulated** (padj < 0.05, log2FC ≤ −1) | **1,738** |

- **Up:** E2F targets (FDR 3×10⁻⁸¹), G2M checkpoint, MYC targets, DNA replication and repair. Top genes include *DLL3, CKAP2L, NDC80, KIF20A*.
- **Down:** phototransduction and sensory perception. Photoreceptor genes such as *NRL, GRK1, RHO* and *GNAT1* are lower.
- **Sanity check:** *RB1* down (log2FC −1.68), *MYCN* up (+4.55), *MKI67* and *TOP2A* up, as expected for retinoblastoma.
- **Interpretation:** without RB1, E2F drives the cell cycle, so proliferation genes switch on while specialised retinal genes drop.

## Figures
| Figure | Description |
|---|---|
| [`volcano_plot.png`](figures/volcano_plot.png) | Volcano plot with gene symbols |
| [`enrichment_heatmap.png`](figures/enrichment_heatmap.png) | Hallmark, KEGG and Reactome enrichment for the up and down gene lists |
| [`enrichment_dotplots.png`](figures/enrichment_dotplots.png) | Top GO Biological Process and Reactome terms |
| [`enrichment_plots.png`](figures/enrichment_plots.png) | Enrichment summary (over-representation and GSEA) |
| [`gsea_running_scores.png`](figures/gsea_running_scores.png) | GSEA running-score curves |
| [`pathway_gene_heatmap.png`](figures/pathway_gene_heatmap.png) | log2FC of the strongest genes in key pathways |
| [`MA_plot.png`](figures/MA_plot.png) | Mean expression versus fold change |
| [`top_genes_and_qc.png`](figures/top_genes_and_qc.png) | Top genes, p-value histogram, gene counts |
| [`pca_samples.png`](figures/pca_samples.png) | PCA of the 7 samples in the GEO raw counts file (74.6% of variance on PC1) |
| [`heatmap_key_genes.png`](figures/heatmap_key_genes.png) | log2 CPM of proliferation, RB1/MYCN and photoreceptor genes per sample |
| [`heatmap_top_DE_genes_samples.png`](figures/heatmap_top_DE_genes_samples.png) | Top 20 up and top 20 down DE genes across the 7 samples (row z-score) |

**Interactive volcano plot:** [`docs/volcano_interactive.html`](docs/volcano_interactive.html). GitHub shows HTML files as code, so to view it, turn on GitHub Pages (Settings → Pages → Branch `main`, folder `/docs`) and open `https://<your-username>.github.io/<repo-name>/volcano_interactive.html`.

![Enrichment](figures/enrichment_plots.png)

## SQL
Results are loaded into a SQLite database ([`data/retinoblastoma.db`](data/retinoblastoma.db)) with four tables: `genes`, `enrichment_ora`, `enrichment_gsea`, `hallmark_membership`. Seven queries ([`sql/analysis_queries.sql`](sql/analysis_queries.sql)) use joins, CTEs, window functions, aggregation with `HAVING`, and `CASE` binning. Their outputs are in [`sql/query_results.md`](sql/query_results.md).

Example: which Hallmark pathways have the highest mean fold change?
```sql
SELECT m.pathway,
       COUNT(*)                AS genes_tested,
       SUM(g.status = 'Up')    AS n_up,
       ROUND(AVG(g.log2fc), 2) AS mean_log2fc
FROM hallmark_membership m
JOIN genes g ON g.symbol = m.symbol
GROUP BY m.pathway
HAVING COUNT(*) >= 30
ORDER BY mean_log2fc DESC
LIMIT 5;
```
| pathway | genes_tested | n_up | mean_log2fc |
|---|---|---|---|
| HALLMARK_E2F_TARGETS | 198 | 135 | 1.99 |
| HALLMARK_G2M_CHECKPOINT | 198 | 116 | 1.77 |
| HALLMARK_MYC_TARGETS_V2 | 56 | 16 | 0.94 |
| HALLMARK_MYC_TARGETS_V1 | 199 | 60 | 0.80 |
| HALLMARK_SPERMATOGENESIS | 101 | 26 | 0.73 |

**To run the queries without coding:** install the free [DB Browser for SQLite](https://sqlitebrowser.org/), open `data/retinoblastoma.db`, go to *Execute SQL*, and paste a query from `sql/analysis_queries.sql`.

## Methods
| Step | Detail |
|---|---|
| Filtering | Genes with a valid adjusted p-value (zero-count, outlier and low-expression genes are removed by DESeq2) |
| Significance | padj < 0.05 and \|log2FC\| ≥ 1 |
| Annotation | Ensembl ID to gene symbol via the MSigDB v7.5.1 mapping (89% of tested genes mapped) |
| Over-representation | One-sided hypergeometric test on the up and down lists separately, background = tested genes, gene sets of 10 to 500 genes, Benjamini-Hochberg FDR |
| GSEA | Preranked on the Wald statistic, 1,000 permutations, `gseapy` |
| Gene sets | MSigDB Hallmark, KEGG, Reactome, GO Biological Process |

## Repository structure
```
├── README.md
├── REPORT.pdf                full written report
├── data/
│   ├── retinoblastoma_gene_results.csv   input DESeq2 table
│   ├── annotated_results_all_genes.csv   input + gene symbols + category
│   ├── significant_DEGs_annotated.csv    3,906 significant genes
│   ├── enrichment_ORA.csv, enrichment_GSEA.csv
│   ├── top15_up.csv, top15_down.csv
│   └── retinoblastoma.db                 SQLite database
├── figures/                  all figures (PNG/PDF)
├── sql/                      queries and their outputs
├── scripts/                  analysis code (see scripts/README.md)
└── docs/                     interactive plot and LinkedIn slides
```

## Reproducing
```bash
pip install -r requirements.txt
python scripts/reproduce_volcano.py              # re-creates the volcano plot from the input table
python scripts/pca_and_sample_heatmaps.py        # PCA and sample heatmaps from the GEO raw counts
```
The other scripts are included for transparency. They were run in a different environment, and the enrichment scripts need the MSigDB gene-set tables, which are not redistributed here. See [`scripts/README.md`](scripts/README.md).

## Limitations
- **Tissue composition (main caveat):** the comparison group is adult retina, which is rich in photoreceptors, while the tumour is a proliferating mass. The drop in photoreceptor and phototransduction genes therefore largely reflects tissue composition and should not be read as tumour cells actively switching those genes off.
- **Contrast direction** is not stored in the input table. Positive log2FC is taken to mean higher in tumour, which is consistent with *RB1* down and *MYCN* up. Confirm against the original GEO file.
- **Annotation:** 11% of tested genes had no symbol and were left out of enrichment.
- **Noisy fold changes:** log2FC for low-count genes is noisy (median lfcSE about 0.6). Applying `lfcShrink` would help before ranking by effect size.
- **Sample labels:** the GEO raw counts file (`data/GSE125903_raw_counts_GRCh38_p13_NCBI.tsv.gz`) has 7 samples and no tumour/retina labels, while the GEO series has 10. Two samples (GSM3584494, GSM3584495) show a retina-like, photoreceptor-high and low-proliferation profile and separate from the other five on PC1, so the sample figures group samples by expression pattern, not by official labels. Check the sample list on the GEO page. Counts for the remaining three samples were not in this file, and the DESeq2 table may have used all 10, so the sample plots do not reproduce the DESeq2 contrast exactly.
- This is a **re-analysis that recovers known biology**, not a new discovery.

## How this was made
The Python analysis scripts and SQL queries were generated with AI assistance (Claude, Anthropic). The project owner framed the question, supplied the DESeq2 results table, reviewed the outputs, checked them against known retinoblastoma biology and wrote the interpretation. The differential expression statistics (DESeq2) were produced before this repository. Please re-run and verify the code before reusing it.

## Data and tools
- Input: DESeq2 results table derived from public dataset **GEO [GSE125903](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125903)** (retinoblastoma tumours vs adult retina; Illumina HiSeq 2500). Original study: Rajasekaran S, Nagarajha Selvan LD, Dotts K, Kumar R, et al. Non-coding and coding transcriptional profiles are significantly altered in pediatric retinoblastoma tumors. *Front Oncol* 2019;9:221 (PMID 31058073). The data belong to the original authors and are shared here with citation for educational and portfolio purposes.
- Tools: Python (pandas, SciPy, matplotlib, seaborn, Plotly, gseapy), SQLite, MSigDB v7.5.1 (Liberzon et al.).

## License
Code: MIT (see `LICENSE`). Check the original data's licence before reuse.
