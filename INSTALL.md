# Hermes — Installation & Database Setup Guide

Complete setup guide for a fresh machine. Follow sections in order.
Estimated total time: 2–4 days (dominated by NR database download ~500 GB).

---

## 1. System Requirements

| Item | Minimum | Recommended |
|---|---|---|
| CPU cores | 16 | 64–160 |
| RAM | 64 GB | 256 GB+ |
| Storage | 2 TB | 4 TB+ (NVMe preferred) |
| GPU | — | NVIDIA RTX with CUDA 12+ (for wgs_snp DeepVariant) |
| OS | Ubuntu 20.04+ | Ubuntu 24.04 |

---

## 2. Base Software

### 2.1 Miniforge3 (conda/mamba)

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash
source ~/.bashrc
```

### 2.2 Snakemake

```bash
mamba install -c conda-forge -c bioconda snakemake
```

### 2.3 Docker (required for wgs_snp DeepVariant GPU)

```bash
curl -fsSL https://get.docker.com | sudo sh
sudo usermod -aG docker $USER
newgrp docker
# Verify
docker run hello-world
```

### 2.4 NVIDIA Container Toolkit (required for wgs_snp GPU)

```bash
distribution=$(. /etc/os-release; echo $ID$VERSION_ID)
curl -s -L https://nvidia.github.io/libnvidia-container/gpgkey | sudo apt-key add -
curl -s -L https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list \
    | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list
sudo apt-get update && sudo apt-get install -y nvidia-container-toolkit
sudo systemctl restart docker
# Verify
docker run --rm --gpus all nvidia/cuda:12.0-base-ubuntu22.04 nvidia-smi
```

---

## 3. Clone Hermes

```bash
git clone https://github.com/cylin2022/hermes.git
cd hermes
```

---

## 4. Database Downloads

Set the base path first; adjust to your storage location:

```bash
export DB_BASE="/path/to/databases"   # e.g. /data/databases
mkdir -p "$DB_BASE"
```

---

### 4.1 NCBI NR — DIAMOND protein database

Used by: `genome_annotation` (functional annotation)
Size: ~300 GB compressed → ~500 GB uncompressed + ~150 GB .dmnd index
Time: 1–2 days download + ~12 h diamond makedb (64 threads)

```bash
mkdir -p "$DB_BASE/nr"
cd "$DB_BASE/nr"

# Download NR FASTA (split into ~8 files, download in parallel)
for i in $(seq -w 00 07); do
    wget -c "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz" \
         -O "nr_part${i}.gz" &
done
# Single file (simpler, slower):
# wget -c "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"

wait
cat nr_part*.gz | gunzip > nr.faa   # or: gunzip -c nr.gz > nr.faa

# Build DIAMOND database (use --threads = available cores)
diamond makedb \
    --in nr.faa \
    --db "$DB_BASE/nr.dmnd" \
    --threads 64 \
    --taxonmap /dev/null \
    --log

# Update config_template.yaml:
# nr_diamond_db: "/path/to/databases/nr.dmnd"
```

---

### 4.2 InterProScan — domain / GO / pathway annotation

Used by: `genome_annotation`
Size: ~30 GB
Time: ~2 h download

```bash
IPR_VER="5.73-104.0"
mkdir -p "$DB_BASE"
cd "$DB_BASE"

wget -c "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${IPR_VER}/interproscan-${IPR_VER}-64-bit.tar.gz"
tar -xzf "interproscan-${IPR_VER}-64-bit.tar.gz"

# Test installation
"$DB_BASE/interproscan-${IPR_VER}/interproscan.sh" --version

# Update config_template.yaml:
# interproscan_dir: "/path/to/databases/interproscan-5.73-104.0"
```

---

### 4.3 eggNOG-mapper — COG / KEGG / OG annotation

Used by: `genome_annotation` (optional step)
Size: ~50 GB
Time: ~3 h download

```bash
mamba install -c conda-forge -c bioconda eggnog-mapper

mkdir -p "$DB_BASE/eggnog_data"
download_eggnog_data.py \
    --data_dir "$DB_BASE/eggnog_data" \
    -y

# Update config_template.yaml:
# eggnog_data_dir: "/path/to/databases/eggnog_data"
```

---

### 4.4 FCS-GX — genome contamination screen

Used by: `genome_annotation` (optional pre-annotation QC)
Size: ~100 GB
Time: ~4 h download

```bash
mkdir -p "$DB_BASE/fcs_gx"
cd "$DB_BASE/fcs_gx"

# Download the GXI database files (NCBI FTP)
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.gxi"
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.gxs"

# Pull FCS-GX container
docker pull ncbi/fcs-gx:latest

# Test
docker run --rm ncbi/fcs-gx:latest python3 /app/bin/run_fcsadaptor.py --help
```

---

### 4.5 BUSCO lineage datasets (auto-downloaded on first run)

BUSCO downloads the selected lineage database automatically when first run.
If offline / air-gapped, pre-download manually:

**IMPORTANT — version compatibility:**
- odb10 datasets → requires BUSCO v5 (`busco=5.7.1` in `envs/busco.yaml`)
- odb12 datasets → requires BUSCO v6 (`busco=6.0.0` in `envs/busco.yaml`)

The `envs/busco.yaml` currently uses **BUSCO v6.0.0** (odb12-compatible).

```bash
# Example: vertebrata (odb12, requires BUSCO v6)
busco --download vertebrata_odb12 --download_path "$DB_BASE/busco_lineages"

# Other common lineages (odb12):
# actinopterygii_odb12  (ray-finned fish)
# mollusca_odb12        (mollusks)
# embryophyta_odb12     (plants)
# insecta_odb12
# mammalia_odb12
```

---

### 4.6 Kraken2 Standard Database — taxonomic profiling

Used by: `metagenome` (Module A)
Size: ~85 GB
Time: ~4 h download + ~30 min build

```bash
mkdir -p "$DB_BASE/kraken2_db"

# Option A: pre-built standard database (fastest)
wget -c "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz" \
    -O "$DB_BASE/kraken2_db/k2_standard.tar.gz"
tar -xzf "$DB_BASE/kraken2_db/k2_standard.tar.gz" -C "$DB_BASE/kraken2_db/"

# Option B: build from scratch (ensures latest taxonomy)
# mamba install -c bioconda kraken2 bracken
# kraken2-build --standard --db "$DB_BASE/kraken2_db" --threads 64
# bracken-build -d "$DB_BASE/kraken2_db" -t 64 -l 150

# Verify
ls "$DB_BASE/kraken2_db"   # should show hash.k2d  opts.k2d  taxo.k2d

# Update metagenome config:
# kraken2_db: "/home/cylin/Vet_Hamaguri/databases/kraken2_db"
```

---

### 4.7 CheckM2 Database — MAG quality assessment

Used by: `metagenome` (Module B)
Size: ~3.5 GB
Time: ~20 min download

```bash
mkdir -p "$DB_BASE/checkm2_db"
mamba install -c bioconda checkm2

# Download using the built-in command
checkm2 database --download --path "$DB_BASE/checkm2_db"

# Verify
ls "$DB_BASE/checkm2_db"   # should show uniref100.KO.1.dmnd

# Update metagenome config:
# checkm2_db: "/home/cylin/Vet_Hamaguri/databases/checkm2_db"
```

---

### 4.8 GTDB-Tk r220 Database — MAG taxonomy

Used by: `metagenome` (Module B)
Size: ~110 GB
Time: ~6 h download

```bash
mkdir -p "$DB_BASE/gtdbtk_db"
mamba install -c bioconda gtdbtk

# Download r220 (latest as of 2024)
wget -c "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz" \
    -O "$DB_BASE/gtdbtk_db/gtdbtk_r220_data.tar.gz"
tar -xzf "$DB_BASE/gtdbtk_db/gtdbtk_r220_data.tar.gz" -C "$DB_BASE/gtdbtk_db/" --strip-components 1

# Set env variable for manual runs (Snakemake workflow sets it automatically)
export GTDBTK_DATA_PATH="$DB_BASE/gtdbtk_db"

# Verify
gtdbtk check_install   # should print 'gtdbtk v2.x.x — OK'

# Update metagenome config:
# gtdbtk_db: "/home/cylin/Vet_Hamaguri/databases/gtdbtk_db"
```

---

### 4.9 Human T2T reference (optional — gut metagenome host removal)

Used by: `metagenome` (host_removal step, only for samples with host_genome set)
Size: ~3 GB
Time: ~10 min download

```bash
mkdir -p "$DB_BASE/host_genomes"
wget -c "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz" \
    -O "$DB_BASE/host_genomes/chm13v2.0.fa.gz"
gunzip "$DB_BASE/host_genomes/chm13v2.0.fa.gz"

# Add to samples.csv for human gut samples:
# gut_01,/path/to/gut_01.hifi.fastq.gz,/home/cylin/Vet_Hamaguri/databases/host_genomes/chm13v2.0.fa
# Leave host_genome blank for soil/environmental samples
```

---

## 5. DeepVariant Docker Image (wgs_snp)

Use the GPU image (`1.10.0-gpu` or later). The older `1.6.1` image exits with code 1 even on success, causing Snakemake to delete all outputs in an infinite failure loop. `1.10.0-gpu` also internalizes tabix indexing, so the Snakefile uses `tabix -f` to force-overwrite the index.

```bash
docker pull google/deepvariant:1.10.0-gpu
# Verify GPU access inside container
docker run --rm --gpus all google/deepvariant:1.10.0-gpu \
    python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

---

## 6. Platform Auto-Configuration (run once per machine)

After all software and databases are installed, run the platform configurator to detect your hardware and generate an optimized config for this machine:

```bash
cd /path/to/hermes

# Quick check (skips GPU container test, takes ~2s)
python3 hermes_configure.py --skip-gpu-test

# Full check including GPU passthrough validation (~30s)
python3 hermes_configure.py

# Write optimized config directly to a file
python3 hermes_configure.py --workflow wgs_snp --output my_server_config.yaml
```

The script reports:
- CPU cores, RAM, GPU model/VRAM
- Whether Docker is installed and the current user is in the `docker` group
- Whether the NVIDIA Container Toolkit can pass GPUs into containers
- Recommended `threads`, `dv_shards`, and `mem_gb` values for this machine

**Key insight**: `threads` controls CPU-based tools (BWA-MEM2, samtools). `dv_shards` controls DeepVariant's internal parallelism. They are set independently so GPU jobs never wait for CPU slots to be free.

If Docker group membership fails:
```bash
sudo usermod -aG docker $USER
newgrp docker        # or log out and back in
```

## 7. MCP Server Setup (Claude Code integration)

The MCP server allows Claude Code to launch and monitor workflows via conversation.

```bash
# Install Python dependencies for the MCP server
pip install mcp snakemake

# Register with Claude Code (add to ~/.claude/settings.json or project settings):
# {
#   "mcpServers": {
#     "hermes": {
#       "command": "python3",
#       "args": ["/path/to/hermes/mcp_server.py"]
#     }
#   }
# }
```

Restart Claude Code after editing settings. Verify by asking Claude:
> "list available workflows"

---

## 7. Summary: config paths by workflow

After all downloads, update each workflow's `config_template.yaml` with the actual paths:

| Config key | Set in | Points to |
|---|---|---|
| `nr_diamond_db` | `genome_annotation/config_template.yaml` | `$DB_BASE/nr.dmnd` |
| `interproscan_dir` | `genome_annotation/config_template.yaml` | `$DB_BASE/interproscan-5.73-104.0` |
| `eggnog_data_dir` | `genome_annotation/config_template.yaml` | `$DB_BASE/eggnog_data` |
| `dv_docker_image` | `wgs_snp/config_template.yaml` | `google/deepvariant:1.10.0-gpu` |
| `kraken2_db` | `metagenome/config.yaml` | `$DB_BASE/kraken2_db` |
| `checkm2_db` | `metagenome/config.yaml` | `$DB_BASE/checkm2_db` |
| `gtdbtk_db` | `metagenome/config.yaml` | `$DB_BASE/gtdbtk_db` |

---

## 8. First run test (rnaseq — no large databases needed)

```bash
cd hermes

# Dry-run with the template (will fail on missing paths — expected)
snakemake -s workflows/rnaseq/Snakefile \
    --use-conda \
    --configfile workflows/rnaseq/config_template.yaml \
    --cores 8 \
    --dry-run
```

If Snakemake prints a job list without errors, the installation is correct.

---

## 9. Storage budget summary

| Database | Compressed | On-disk |
|---|---|---|
| NR FASTA + DIAMOND index | ~300 GB + ~150 GB | ~450 GB |
| InterProScan | — | ~30 GB |
| eggNOG | — | ~50 GB |
| FCS-GX | — | ~100 GB |
| BUSCO lineages (all) | — | ~5 GB |
| DeepVariant image | — | ~3 GB |
| Kraken2 Standard DB | — | ~85 GB |
| CheckM2 DB | — | ~3.5 GB |
| GTDB-Tk r220 | — | ~110 GB |
| Human T2T (optional) | — | ~3 GB |
| **Total (all)** | | **~840 GB** |
