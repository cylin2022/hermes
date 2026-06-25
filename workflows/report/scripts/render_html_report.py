"""
Universal HTML report renderer.
Reads all input files, embeds tables and base64-encoded figures into a
self-contained HTML file. The template is selected by snakemake.config["workflow_type"].
"""
import base64, os, re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from jinja2 import Template

WORKFLOW = snakemake.config["workflow_type"]
PROJECT  = snakemake.params.project
SAMPLES  = snakemake.params.samples
OUT      = snakemake.output[0]
INPUT    = dict(snakemake.input)   # named inputs as dict

def img_b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()

def df_to_html(path, max_rows=200):
    df = pd.read_csv(path, sep="\t")
    return df.head(max_rows).to_html(
        index=False, classes="tbl", border=0,
        float_format=lambda x: f"{x:.3f}",
    )

# ── Collect section data by workflow type ──────────────────────────────────────
sections = []

if WORKFLOW == "metagenome":
    if "mag_summary" in INPUT:
        df = pd.read_csv(INPUT["mag_summary"], sep="\t")
        tiers = df["quality_tier"].value_counts().to_dict() if "quality_tier" in df.columns else {}
        sections.append({
            "title": "MAG Quality Summary",
            "text": f"Total MAGs: {len(df)} | " +
                    " | ".join(f"{k}: {v}" for k, v in tiers.items()),
            "table_html": df_to_html(INPUT["mag_summary"]),
            "figure_b64": img_b64(INPUT["scatter"]) if "scatter" in INPUT else None,
            "figure_caption": "Completeness vs. Contamination (MIMAG thresholds shown as dashed lines)",
        })
    if "taxonomy" in INPUT:
        sections.append({
            "title": "Community Taxonomy (Species Level)",
            "text":  "Bracken species-level relative abundance across samples.",
            "table_html": df_to_html(INPUT["taxonomy"]),
            "figure_b64": img_b64(INPUT["barplot"]) if "barplot" in INPUT else None,
            "figure_caption": f"Top species by mean abundance",
        })

elif WORKFLOW == "wgs_snp":
    if "variants" in INPUT:
        sections.append({
            "title": "Variant Call Summary",
            "text":  "SNP/INDEL counts and Ti/Tv ratios.",
            "table_html": df_to_html(INPUT["variants"]),
            "figure_b64": img_b64(INPUT["titv"]) if "titv" in INPUT else None,
            "figure_caption": "Ti/Tv ratio per sample",
        })

elif WORKFLOW == "rnaseq":
    if "de_tbl" in INPUT:
        df = pd.read_csv(INPUT["de_tbl"], sep="\t")
        sig = (df["padj"] < 0.05).sum() if "padj" in df.columns else "N/A"
        sections.append({
            "title": "Differential Expression Summary",
            "text":  f"Significant DEGs (padj < 0.05): {sig}",
            "table_html": df_to_html(INPUT["de_tbl"]),
            "figure_b64": img_b64(INPUT["pca"]) if "pca" in INPUT else None,
            "figure_caption": "PCA of normalised counts",
        })

elif WORKFLOW == "scrnaseq":
    if "clusters" in INPUT:
        sections.append({
            "title": "Cluster Summary",
            "text":  "Leiden clustering marker genes.",
            "table_html": df_to_html(INPUT["clusters"]),
            "figure_b64": img_b64(INPUT["umap"]) if "umap" in INPUT else None,
            "figure_caption": "UMAP coloured by Leiden cluster",
        })

# ── HTML template ──────────────────────────────────────────────────────────────
TMPL = Template(r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ project }} — Hermes Report</title>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         margin: 0; padding: 0; background: #f8f9fa; color: #212529; }
  header { background: #1a1a2e; color: #fff; padding: 24px 40px; }
  header h1 { margin: 0 0 4px; font-size: 1.6em; }
  header p  { margin: 0; opacity: 0.75; font-size: 0.9em; }
  main { max-width: 1200px; margin: 32px auto; padding: 0 24px; }
  section { background: #fff; border-radius: 8px; box-shadow: 0 1px 4px rgba(0,0,0,.08);
            margin-bottom: 32px; padding: 24px 28px; }
  h2 { margin-top: 0; border-bottom: 2px solid #e9ecef; padding-bottom: 8px;
       font-size: 1.2em; color: #343a40; }
  .meta { color: #6c757d; font-size: 0.88em; margin-bottom: 12px; }
  .tbl-wrap { overflow-x: auto; margin: 12px 0; }
  table.tbl { border-collapse: collapse; width: 100%; font-size: 0.82em; }
  table.tbl th { background: #343a40; color: #fff; padding: 6px 10px; text-align: left; }
  table.tbl td { padding: 5px 10px; border-bottom: 1px solid #dee2e6; }
  table.tbl tr:hover td { background: #f1f3f5; }
  figure { margin: 0; text-align: center; }
  figure img { max-width: 100%; border-radius: 6px; box-shadow: 0 1px 6px rgba(0,0,0,.12); }
  figcaption { color: #6c757d; font-size: 0.8em; margin-top: 6px; }
  .multiqc-link { display: inline-block; margin-top: 8px; padding: 8px 16px;
                  background: #0d6efd; color: #fff; border-radius: 4px;
                  text-decoration: none; font-size: 0.9em; }
  .multiqc-link:hover { background: #0b5ed7; }
  footer { text-align: center; color: #adb5bd; font-size: 0.8em; padding: 24px; }
</style>
</head>
<body>
<header>
  <h1>{{ project }}</h1>
  <p>Workflow: <strong>{{ workflow }}</strong> &nbsp;|&nbsp;
     Samples: <strong>{{ n_samples }}</strong> &nbsp;|&nbsp;
     Generated by <strong>小賀 / Hermes</strong></p>
</header>
<main>

  <section>
    <h2>QC Summary</h2>
    <p class="meta">MultiQC aggregates fastp read quality and adapter statistics.</p>
    <a class="multiqc-link" href="{{ multiqc_path }}" target="_blank">Open MultiQC Report ↗</a>
  </section>

  {% for sec in sections %}
  <section>
    <h2>{{ sec.title }}</h2>
    {% if sec.text %}<p class="meta">{{ sec.text }}</p>{% endif %}
    {% if sec.figure_b64 %}
    <figure>
      <img src="data:image/png;base64,{{ sec.figure_b64 }}" alt="{{ sec.figure_caption }}">
      <figcaption>{{ sec.figure_caption }}</figcaption>
    </figure>
    {% endif %}
    {% if sec.table_html %}
    <div class="tbl-wrap">{{ sec.table_html }}</div>
    {% endif %}
  </section>
  {% endfor %}

</main>
<footer>Generated by <a href="https://github.com/cylin2022/hermes">Hermes</a></footer>
</body>
</html>
""")

html = TMPL.render(
    project      = PROJECT,
    workflow     = WORKFLOW,
    n_samples    = len(SAMPLES),
    multiqc_path = INPUT.get("multiqc", "#"),
    sections     = sections,
)

Path(OUT).write_text(html, encoding="utf-8")
print(f"Report written: {OUT}")
