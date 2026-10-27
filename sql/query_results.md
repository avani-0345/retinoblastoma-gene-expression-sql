# SQL query results

Output of each query in `analysis_queries.sql`, run on `data/retinoblastoma.db` (SQLite).


## Q1. Overview: genes per category, with average expression and effect size

```sql
-- Q1. Overview: genes per category, with average expression and effect size
SELECT status,
       COUNT(*)                         AS n_genes,
       ROUND(AVG(base_mean), 1)         AS avg_base_mean,
       ROUND(AVG(ABS(log2fc)), 2)       AS avg_abs_log2fc,
       ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 1) AS pct_of_tested
FROM genes
GROUP BY status
ORDER BY n_genes DESC;
```

| status          |   n_genes |   avg_base_mean |   avg_abs_log2fc |   pct_of_tested |
|:----------------|----------:|----------------:|-----------------:|----------------:|
| NS              |     15793 |           741.7 |             0.77 |            79.2 |
| Up              |      2168 |           996   |             2.84 |            10.9 |
| Down            |      1738 |          1390.6 |             2.19 |             8.7 |
| Sig. (small FC) |       253 |          1901.6 |             0.87 |             1.3 |



## Q2. Top 5 up- and down-regulated genes by adjusted p-value (window function)

```sql
-- Q2. Top 5 up- and down-regulated genes by adjusted p-value (window function)
WITH ranked AS (
    SELECT symbol, ensembl_id, status, ROUND(log2fc, 2) AS log2fc, padj,
           ROW_NUMBER() OVER (PARTITION BY status ORDER BY padj ASC) AS rn
    FROM genes
    WHERE status IN ('Up', 'Down')
)
SELECT status, rn AS rank, symbol, log2fc, padj
FROM ranked
WHERE rn <= 5
ORDER BY status DESC, rn;
```

| status   |   rank | symbol   |   log2fc |     padj |
|:---------|-------:|:---------|---------:|---------:|
| Up       |      1 | DLL3     |     6.31 | 1.68e-38 |
| Up       |      2 | CKAP2L   |     9.23 | 8.8e-38  |
| Up       |      3 | NDC80    |     9.01 | 1.49e-34 |
| Up       |      4 | KIF20A   |     8.28 | 5.24e-33 |
| Up       |      5 | SLC7A3   |     6.34 | 8.06e-33 |
| Down     |      1 | RPS27P25 |    -3.7  | 1.99e-18 |
| Down     |      2 | NRL      |    -3.69 | 8.69e-16 |
| Down     |      3 | TAFA3    |    -3.41 | 1.1e-15  |
| Down     |      4 | MIR22HG  |    -3.19 | 7.04e-14 |
| Down     |      5 | CADM2    |    -3.57 | 2.74e-13 |



## Q3. Known retinoblastoma / retinal genes (CTE + CASE)

```sql
-- Q3. Known retinoblastoma / retinal genes (CTE + CASE)
WITH panel(symbol, role) AS (
    VALUES ('RB1','Tumour suppressor'), ('MYCN','Oncogene'), ('MKI67','Proliferation'),
           ('TOP2A','Proliferation'), ('NRL','Rod photoreceptor'), ('GRK1','Phototransduction'),
           ('RHO','Phototransduction'), ('DLL3','Notch ligand')
)
SELECT p.symbol, p.role, ROUND(g.log2fc, 2) AS log2fc, g.padj,
       CASE WHEN g.log2fc > 0 THEN 'Up' ELSE 'Down' END AS direction
FROM panel p
JOIN genes g ON g.symbol = p.symbol
ORDER BY g.log2fc DESC;
```

| symbol   | role              |   log2fc |     padj | direction   |
|:---------|:------------------|---------:|---------:|:------------|
| DLL3     | Notch ligand      |     6.31 | 1.68e-38 | Up          |
| MYCN     | Oncogene          |     4.55 | 1.07e-05 | Up          |
| TOP2A    | Proliferation     |     3.93 | 0.000854 | Up          |
| MKI67    | Proliferation     |     2.77 | 0.0314   | Up          |
| RB1      | Tumour suppressor |    -1.68 | 0.00872  | Down        |
| GRK1     | Phototransduction |    -3.54 | 4.7e-13  | Down        |
| NRL      | Rod photoreceptor |    -3.69 | 8.69e-16 | Down        |
| RHO      | Phototransduction |    -4.65 | 0.000736 | Down        |



## Q4. Hallmark pathway activity: JOIN genes to pathway membership, then aggregate

```sql
-- Q4. Hallmark pathway activity: JOIN genes to pathway membership, then aggregate
SELECT m.pathway,
       COUNT(*)                                            AS genes_tested,
       SUM(g.status = 'Up')                                AS n_up,
       SUM(g.status = 'Down')                              AS n_down,
       ROUND(AVG(g.log2fc), 2)                             AS mean_log2fc,
       ROUND(100.0 * SUM(g.status = 'Up') / COUNT(*), 1)   AS pct_up
FROM hallmark_membership m
JOIN genes g ON g.symbol = m.symbol
GROUP BY m.pathway
HAVING COUNT(*) >= 30
ORDER BY mean_log2fc DESC
LIMIT 10;
```

| pathway                          |   genes_tested |   n_up |   n_down |   mean_log2fc |   pct_up |
|:---------------------------------|---------------:|-------:|---------:|--------------:|---------:|
| HALLMARK_E2F_TARGETS             |            198 |    135 |        2 |          1.99 |     68.2 |
| HALLMARK_G2M_CHECKPOINT          |            198 |    116 |        6 |          1.77 |     58.6 |
| HALLMARK_MYC_TARGETS_V2          |             56 |     16 |        0 |          0.94 |     28.6 |
| HALLMARK_MYC_TARGETS_V1          |            199 |     60 |        1 |          0.8  |     30.2 |
| HALLMARK_SPERMATOGENESIS         |            101 |     26 |        9 |          0.73 |     25.7 |
| HALLMARK_MITOTIC_SPINDLE         |            195 |     46 |       14 |          0.59 |     23.6 |
| HALLMARK_ALLOGRAFT_REJECTION     |            156 |     21 |       11 |          0.38 |     13.5 |
| HALLMARK_DNA_REPAIR              |            145 |     29 |        6 |          0.35 |     20   |
| HALLMARK_IL6_JAK_STAT3_SIGNALING |             73 |      7 |        2 |          0.28 |      9.6 |
| HALLMARK_MTORC1_SIGNALING        |            195 |     29 |       22 |          0.17 |     14.9 |



## Q5. Top 3 enriched terms per collection and direction (window function over enrichment table)

```sql
-- Q5. Top 3 enriched terms per collection and direction (window function over enrichment table)
SELECT collection, direction, term, overlap, set_size, fdr
FROM (
    SELECT *, ROW_NUMBER() OVER (PARTITION BY collection, direction ORDER BY fdr ASC) AS rn
    FROM enrichment_ora
    WHERE fdr < 0.05
) t
WHERE rn <= 3
ORDER BY collection, direction DESC, rn;
```

| collection   | direction   | term                                      |   overlap |   set_size |      fdr |
|:-------------|:------------|:------------------------------------------|----------:|-----------:|---------:|
| GO_BP        | Up          | GOBP_CHROMOSOME_SEGREGATION               |       119 |        320 | 2.37e-32 |
| GO_BP        | Up          | GOBP_NUCLEAR_CHROMOSOME_SEGREGATION       |       103 |        261 | 1.51e-30 |
| GO_BP        | Up          | GOBP_DNA_REPAIR                           |       144 |        495 | 1.14e-26 |
| GO_BP        | Down        | GOBP_SENSORY_PERCEPTION_OF_LIGHT_STIMULUS |        71 |        185 | 1.55e-23 |
| GO_BP        | Down        | GOBP_SENSORY_PERCEPTION                   |       102 |        423 | 7.68e-17 |
| GO_BP        | Down        | GOBP_DETECTION_OF_LIGHT_STIMULUS          |        27 |         58 | 1.55e-10 |
| Hallmark     | Up          | HALLMARK_E2F_TARGETS                      |       135 |        198 | 3.04e-81 |
| Hallmark     | Up          | HALLMARK_G2M_CHECKPOINT                   |       116 |        198 | 3.71e-59 |
| Hallmark     | Up          | HALLMARK_MYC_TARGETS_V1                   |        60 |        199 | 1.32e-12 |
| Hallmark     | Down        | HALLMARK_MYOGENESIS                       |        36 |        153 | 5.29e-06 |
| Hallmark     | Down        | HALLMARK_ESTROGEN_RESPONSE_EARLY          |        36 |        171 | 4.64e-05 |
| Hallmark     | Down        | HALLMARK_XENOBIOTIC_METABOLISM            |        33 |        165 | 0.000252 |
| KEGG         | Up          | KEGG_CELL_CYCLE                           |        49 |        119 | 1.48e-15 |
| KEGG         | Up          | KEGG_DNA_REPLICATION                      |        25 |         35 | 2.57e-15 |
| KEGG         | Up          | KEGG_P53_SIGNALING_PATHWAY                |        24 |         65 | 1.19e-06 |
| KEGG         | Down        | KEGG_ARGININE_AND_PROLINE_METABOLISM      |        14 |         51 | 0.0197   |
| KEGG         | Down        | KEGG_DILATED_CARDIOMYOPATHY               |        17 |         75 | 0.0227   |
| KEGG         | Down        | KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM      |        16 |         70 | 0.0227   |
| Reactome     | Up          | REACTOME_CELL_CYCLE_CHECKPOINTS           |        95 |        263 | 8.11e-25 |
| Reactome     | Up          | REACTOME_MITOTIC_PROMETAPHASE             |        75 |        196 | 2.44e-21 |
| Reactome     | Up          | REACTOME_M_PHASE                          |       104 |        364 | 1.15e-18 |
| Reactome     | Down        | REACTOME_SENSORY_PERCEPTION               |        55 |        177 | 1.32e-13 |
| Reactome     | Down        | REACTOME_VISUAL_PHOTOTRANSDUCTION         |        32 |         81 | 7.73e-11 |
| Reactome     | Down        | REACTOME_THE_PHOTOTRANSDUCTION_CASCADE    |        18 |         29 | 6.04e-10 |



## Q6. Agreement between ORA and GSEA for Hallmark (self-join across tables)

```sql
-- Q6. Agreement between ORA and GSEA for Hallmark (self-join across tables)
SELECT o.term,
       ROUND(g.nes, 2) AS gsea_nes,
       o.direction     AS ora_direction,
       o.fdr           AS ora_fdr,
       g.fdr           AS gsea_fdr
FROM enrichment_ora o
JOIN enrichment_gsea g ON g.term = o.term AND g.collection = 'Hallmark'
WHERE o.collection = 'Hallmark' AND o.fdr < 0.05 AND g.fdr < 0.05
ORDER BY o.fdr
LIMIT 8;
```

| term                             |   gsea_nes | ora_direction   |   ora_fdr |   gsea_fdr |
|:---------------------------------|-----------:|:----------------|----------:|-----------:|
| HALLMARK_E2F_TARGETS             |       3.54 | Up              |  3.04e-81 |    0       |
| HALLMARK_G2M_CHECKPOINT          |       3.35 | Up              |  3.71e-59 |    0       |
| HALLMARK_MYC_TARGETS_V1          |       2.75 | Up              |  1.32e-12 |    0       |
| HALLMARK_MITOTIC_SPINDLE         |       2.05 | Up              |  3.54e-06 |    0       |
| HALLMARK_MYOGENESIS              |      -1.62 | Down            |  5.29e-06 |    0.00799 |
| HALLMARK_ESTROGEN_RESPONSE_EARLY |      -1.68 | Down            |  4.64e-05 |    0.00758 |
| HALLMARK_SPERMATOGENESIS         |       2.08 | Up              |  0.000206 |    0       |
| HALLMARK_XENOBIOTIC_METABOLISM   |      -1.69 | Down            |  0.000252 |    0.00799 |



## Q7. Effect-size distribution of significant genes (binning with CASE)

```sql
-- Q7. Effect-size distribution of significant genes (binning with CASE)
SELECT CASE WHEN ABS(log2fc) < 2 THEN '1-2'
            WHEN ABS(log2fc) < 4 THEN '2-4'
            WHEN ABS(log2fc) < 6 THEN '4-6'
            ELSE '6+' END AS abs_log2fc_bin,
       SUM(status = 'Up')   AS up,
       SUM(status = 'Down') AS down
FROM genes
WHERE status IN ('Up','Down')
GROUP BY abs_log2fc_bin
ORDER BY abs_log2fc_bin;
```

| abs_log2fc_bin   |   up |   down |
|:-----------------|-----:|-------:|
| 1-2              |  799 |    843 |
| 2-4              |  930 |    839 |
| 4-6              |  336 |     56 |
| 6+               |  103 |      0 |

