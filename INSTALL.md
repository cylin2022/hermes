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

```bash
# Example: vertebrata
busco --download vertebrata_odb10 --download_path "$DB_BASE/busco_lineages"

# Other common lineages:
# actinopterygii_odb10  (ray-finned fish)
# embryophyta_odb10     (plants)
# insecta_odb10
# mammalia_odb10
```

---

## 5. DeepVariant Docker Image (wgs_snp)

```bash
docker pull google/deepvariant:1.6.1
# Verify GPU access inside container
docker run --rm --gpus all google/deepvariant:1.6.1 \
    python3 -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

---

## 6. MCP Server Setup (Claude Code integration)

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
| `dv_docker_image` | `wgs_snp/config_template.yaml` | `google/deepvariant:1.6.1` |

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
| **Total** | | **~640 GB** |
