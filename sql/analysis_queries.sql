-- Retinoblastoma DESeq2 results: SQL analysis (SQLite; works in PostgreSQL with minor changes)
-- Tables: genes, enrichment_ora, enrichment_gsea, hallmark_membership

-- Q1. Overview: genes per category, with average expression and effect size
SELECT status,
       COUNT(*)                         AS n_genes,
       ROUND(AVG(base_mean), 1)         AS avg_base_mean,
       ROUND(AVG(ABS(log2fc)), 2)       AS avg_abs_log2fc,
       ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 1) AS pct_of_tested
FROM genes
GROUP BY status
ORDER BY n_genes DESC;

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

-- Q5. Top 3 enriched terms per collection and direction (window function over enrichment table)
SELECT collection, direction, term, overlap, set_size, fdr
FROM (
    SELECT *, ROW_NUMBER() OVER (PARTITION BY collection, direction ORDER BY fdr ASC) AS rn
    FROM enrichment_ora
    WHERE fdr < 0.05
) t
WHERE rn <= 3
ORDER BY collection, direction DESC, rn;

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
