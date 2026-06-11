# RNA-seq Differential Expression Workflow

Bulk RNA-seq analysis: QC → alignment → quantification → DE → pathway

## Steps
1. FastQC + MultiQC (quality control)
2. Trimmomatic (adapter trimming)
3. STAR (genome alignment, 2-pass)
4. featureCounts (gene-level quantification)
5. DESeq2 (differential expression, R)
6. clusterProfiler (GO + KEGG pathway enrichment)

## Supported species
- Mus musculus (GRCm39)
- Homo sapiens (GRCh38)
- Custom (provide genome + GTF)

## Estimated runtime (160 cores, 12 samples)
~4–6 hours
