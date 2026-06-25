/**
 * Hermes Workflow Critique
 *
 * Spawns independent review agents for each specified workflow,
 * then synthesizes findings into a prioritized fix list.
 *
 * Usage:
 *   Workflow({ name: "critique" })                          — review ALL workflows
 *   Workflow({ name: "critique", args: ["pool_seq"] })      — review one workflow
 *   Workflow({ name: "critique", args: ["pool_seq","snp_association"] }) — review specific list
 *
 * Each agent reads the Snakefile, scripts, and conda envs independently
 * (no shared context) and checks against PIPELINE_CHECKLIST.md criteria.
 */

export const meta = {
  name: 'critique',
  description: 'Independent code review of Hermes bioinformatics workflows before production runs',
  whenToUse: 'Run before committing any new or modified workflow. Pass workflow name(s) as args, or omit to review all.',
  phases: [
    { title: 'Critique', detail: 'Parallel independent review of each workflow' },
    { title: 'Synthesize', detail: 'Prioritized fix list with verdict per workflow' },
  ],
}

// ── Discover which workflows to review ────────────────────────────────────────
const ALL_WORKFLOWS = [
  'pool_seq', 'snp_association', 'genomic_prediction',
  'atacseq', 'scrnaseq', 'metagenome', 'spatial', 'report',
  'rnaseq', 'wgs_snp', 'genome_annotation',
]

const targets = (Array.isArray(args) && args.length > 0)
  ? args
  : ALL_WORKFLOWS

log(`Reviewing ${targets.length} workflow(s): ${targets.join(', ')}`)

// ── Structured output schema ───────────────────────────────────────────────────
const ISSUE_SCHEMA = {
  type: 'object',
  properties: {
    workflow:        { type: 'string' },
    critical_issues: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          category:      { type: 'string', description: 'docker|shell_escape|samtools_sort|io_format|missing_tool|stats|snakemake_rule|conda_env|other' },
          location:      { type: 'string', description: 'filename:line or rule name' },
          description:   { type: 'string' },
          suggested_fix: { type: 'string' },
        },
        required: ['category', 'location', 'description', 'suggested_fix'],
      },
    },
    minor_issues: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          category:      { type: 'string' },
          location:      { type: 'string' },
          description:   { type: 'string' },
          suggested_fix: { type: 'string' },
        },
        required: ['category', 'location', 'description', 'suggested_fix'],
      },
    },
    verdict: { type: 'string', enum: ['ready_to_run', 'minor_fixes_needed', 'critical_fixes_required'] },
    summary: { type: 'string' },
  },
  required: ['workflow', 'critical_issues', 'minor_issues', 'verdict', 'summary'],
}

// ── Critique prompt (one per workflow) ────────────────────────────────────────
const critiquePrompt = (name) => `
You are a senior bioinformatics engineer doing an independent code review of the
"${name}" Snakemake workflow at /home/cylin/hermes/workflows/${name}/.

You have NO knowledge of why design decisions were made. Review with fresh eyes.
Your job: find bugs that would cause crashes or wrong results in production.
The cost of a missed bug is days of lost compute time on a 160-core / 2.2 TiB server.

## Read these files first
- /home/cylin/hermes/workflows/${name}/Snakefile
- /home/cylin/hermes/workflows/${name}/config_template.yaml
- All /home/cylin/hermes/workflows/${name}/envs/*.yaml
- All /home/cylin/hermes/workflows/${name}/scripts/*.R and *.py
- /home/cylin/hermes/PIPELINE_CHECKLIST.md  (the project checklist — check every gate)

Reference implementations (production-tested, use as comparison):
- /home/cylin/hermes/workflows/wgs_snp/Snakefile    (Docker + samtools patterns)
- /home/cylin/hermes/workflows/genome_annotation/Snakefile  (conda + Docker mix)

## What to check (from PIPELINE_CHECKLIST.md — all gates)

**Docker rules (CRITICAL if violated):**
- Every docker run must have --user $(id -u):$(id -g)
  EXCEPTION: BRAKER3 must run as root (no --user)
- mkdir -p output dir BEFORE docker run
- Use realpath for bind-mount paths
- busybox chown to fix root-owned outputs inside shell block

**Shell block rules (CRITICAL):**
- set -euo pipefail as FIRST line of every shell block
- test -s {output.X} || { echo "ERROR" >&2; exit 1; } after every main command
- sed back-references: use \\\\1 not \\1 inside Snakemake shell blocks
  (Python processes the string first: \\1 → SOH control char, not \\1 for sed)
- || true only with an explanatory comment

**samtools order (CRITICAL — silent wrong results):**
- If fixmate is used: input MUST be name-sorted (queryname)
  Correct: bwa-mem2 | samtools sort -n | samtools fixmate -m | samtools sort | samtools markdup
  Wrong: bwa-mem2 | samtools fixmate (multi-thread output not guaranteed name-sorted)

**Index files:**
- .fai, .bai/.csi, .tbi must be declared as inputs where tools need them
- samtools >= 1.12 writes .csi (not .bai) for large/non-standard references

**Conda environment completeness:**
- Every library() in R scripts must have a package in the conda yaml
- Every import in Python scripts must have a package in the conda yaml
- Tool version numbers consistent between Snakefile and config_template.yaml

**Statistical correctness:**
- Normalization method appropriate for data type
- Multiple testing correction applied (p.adjust / FDR)
- No sample name / path hard-coding
- Edge-case guards: zero genes/SNPs, empty outputs, NA handling

## Output
Use StructuredOutput. Classify as CRITICAL (causes crash or wrong result)
vs MINOR (style, sub-optimal, won't fail). Include file:line for every issue.
`

// ── Phase 1: parallel critique ─────────────────────────────────────────────────
phase('Critique')

const results = await parallel(targets.map(name => () =>
  agent(critiquePrompt(name), {
    label:  `critique:${name}`,
    phase:  'Critique',
    schema: ISSUE_SCHEMA,
    effort: 'high',
  })
))

const valid = results.filter(Boolean)
const nCritical = valid.reduce((s, r) => s + r.critical_issues.length, 0)
const nMinor    = valid.reduce((s, r) => s + r.minor_issues.length,    0)
log(`Done: ${valid.length}/${targets.length} reviewed — ${nCritical} critical, ${nMinor} minor`)

// ── Phase 2: synthesis ────────────────────────────────────────────────────────
phase('Synthesize')

const details = valid.map(r => `
=== ${r.workflow} [${r.verdict}] ===
CRITICAL (${r.critical_issues.length}):
${r.critical_issues.map(i =>
  `  [${i.category}] ${i.location}\n    ${i.description}\n    FIX: ${i.suggested_fix}`
).join('\n') || '  (none)'}
MINOR (${r.minor_issues.length}):
${r.minor_issues.map(i =>
  `  [${i.category}] ${i.location}: ${i.description}`
).join('\n') || '  (none)'}
Summary: ${r.summary}
`).join('\n')

const synthesis = await agent(`
You are synthesizing code-review results for Hermes bioinformatics workflows.

${details}

Write a concise, actionable report for the developer:

## 1. Overall verdict (1-2 sentences)
## 2. Workflows ready to run
   List workflows with verdict "ready_to_run" or only minor issues.
## 3. Critical fixes required
   Group by workflow. For each: location, problem, exact fix.
   Order by severity (crashes first, wrong results second).
## 4. Systemic patterns
   Issues appearing in 2+ workflows → propose a single bulk fix.
## 5. Next steps
   Concrete next action the developer should take.

Be specific. Give exact code where possible. Skip prose filler.
`, {
  label:  'synthesis',
  phase:  'Synthesize',
  effort: 'high',
})

return { findings: valid, synthesis }
