# Hermes — Bioinformatics Analysis Agent

Autonomous bioinformatics analysis platform running on a high-performance server (160 cores, 2.2 TiB RAM, Ubuntu 24.04).

## Overview

Hermes provides a set of Snakemake workflows orchestrated via an MCP server, allowing Claude to plan, execute, and monitor bioinformatics analyses end-to-end.

## Available Workflows

| Workflow | Analysis Type | Key Tools |
|----------|--------------|-----------|
| `genome_annotation` | Genome QC + Repeat + BRAKER4 + Functional | BUSCO, RepeatModeler, BRAKER4, DIAMOND, InterPro |
| `rnaseq` | RNA-seq differential expression | STAR, featureCounts, DESeq2, clusterProfiler |
| `atacseq` | ATAC-seq peak calling | Bowtie2, MACS3, DiffBind |
| `scrnaseq` | Single-cell RNA-seq clustering + markers | STARsolo, scDblFinder, Scanpy, Leiden |
| `wgs_snp` | WGS SNP/INDEL calling (non-model organisms) | fastp, BWA-MEM2, DeepVariant (GPU), GLnexus, SNPeff |
| `snp_association` | SNP-trait GWAS + Fst scan | PLINK2, GEMMA LMM, VCFtools Fst, CMplot |
| `pool_seq` | Pool-seq allele frequency + Fst scan | fastp, BWA-MEM2, bcftools, Hudson Fst (R) |

## Input Format

All workflows use standardized CSV samplesheets. See `samplesheet_template.csv` in each workflow directory.

### Sequence workflows (`wgs_snp`, `scrnaseq`, `rnaseq`, `atacseq`)
```csv
sample,r1,r2
sample_001,/path/to/sample_001_R1.fastq.gz,/path/to/sample_001_R2.fastq.gz
```

### GWAS metadata (`snp_association`)
```csv
sample_id,phenotype,weight_g,length_cm
fish_001,1,52.3,16.1
fish_031,0,41.7,14.5
```

## Structure

```
hermes/
├── mcp_server.py          # MCP server (tools: run_workflow, get_status, etc.)
├── setup_hermes.sh        # Environment setup script
├── CLAUDE.md              # Agent instructions
└── workflows/
    ├── wgs_snp/
    ├── snp_association/
    ├── rnaseq/
    ├── scrnaseq/
    ├── atacseq/
    └── genome_annotation/
```

## Key Paths (server-specific)

```
Databases     : /home/cylin/Vet_Hamaguri/databases/
FCS-GX DB     : /home/cylin/Vet_Hamaguri/FCS_DB
Workflow runs : /home/cylin/hermes/runs/<run_id>/
```

## Requirements

- Snakemake ≥ 8
- Conda / Mamba
- Docker (for DeepVariant GPU)
- NVIDIA GPU + CUDA (for `wgs_snp` DeepVariant step)
- Claude Code with MCP server configured
