"""
Split a genome FASTA into N roughly equal-size chunks for parallel RepeatModeler.
Keeps contigs intact; uses greedy bin-packing to balance total bp per chunk.
Called as a Snakemake script (snakemake.input, snakemake.output, snakemake.params).
"""

from pathlib import Path


def read_fasta(path):
    seqs = []
    with open(path) as fh:
        header = seq_parts = None
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    seq = "".join(seq_parts)
                    seqs.append((header, seq, len(seq)))
                header, seq_parts = line, []
            else:
                seq_parts.append(line)
        if header:
            seq = "".join(seq_parts)
            seqs.append((header, seq, len(seq)))
    return seqs


def greedy_split(seqs, n):
    """Assign sequences to n bins minimising max-bin size (greedy longest-first)."""
    seqs_sorted = sorted(seqs, key=lambda x: x[2], reverse=True)
    bins = [[] for _ in range(n)]
    sizes = [0] * n
    for seq in seqs_sorted:
        i = sizes.index(min(sizes))
        bins[i].append(seq)
        sizes[i] += seq[2]
    return bins


genome_fa = snakemake.input.genome
n_chunks  = snakemake.params.n
outdir    = Path(snakemake.params.outdir)
outdir.mkdir(parents=True, exist_ok=True)

seqs  = read_fasta(genome_fa)
bins  = greedy_split(seqs, n_chunks)

for i, chunk in enumerate(bins):
    outpath = outdir / f"chunk_{i:03d}.fa"
    with open(outpath, "w") as fh:
        if not chunk:
            fh.write(f">chunk_{i}_placeholder\nN\n")
        else:
            for header, seq, _ in chunk:
                fh.write(f"{header}\n")
                for j in range(0, len(seq), 80):
                    fh.write(seq[j:j+80] + "\n")
