# SV–environment workflow with PC3, PC5 and PC10

This workflow starts from **observed SV genotypes** and a frozen sample-to-environment join. It computes one nested SV PCA and fits LFMM2 and partial RDA at 3, 5 and 10 explicit PC covariates. It uses no SNP input. A fixed LFMM2 latent-factor count is independent of the PC count. See the [Chinese methods and operations manual](docs/SV_PC_LFMM2_RDA_操作与原理.pdf).

## Quick start

Use Linux, Python 3.10+ with NumPy and pandas, and R with LEA, vegan, permute and jsonlite. Use Bioconductor to install LEA. VCF/BCF input requires bcftools. The example is fictional and the complete run typically takes seconds to tens of seconds on a workstation.

```bash
python3 -m pip install -r requirements.txt
python3 scripts/workflow.py example/run.json runs/synthetic_demo
bash status.sh runs/synthetic_demo
```

Set `SV_RSCRIPT` to a custom Rscript executable or `SV_R_LIBRARY` to its additional R library path. For a long run, set `SV_PYTHON` as needed and use `bash start.sh CONFIG.json OUTPUT_DIRECTORY`. The detached worker runs independently of the tmux monitor. Choose a fresh output directory when changing inputs, code, software or parameters.

## Formal inputs

Copy `configs/run.template.json` and update paths, variables and settings. Relative input paths are interpreted relative to the JSON file. Use the verified project cohort metadata to create these tables.

| Input | Contract |
|---|---|
| SV genotype | VCF/BCF with observed GT calls or a sample-by-locus TSV with a first column named `sample` and diploid 0/1/2 or NA |
| Sample table | TSV with unique `sample` and a `group` for each included individual |
| Group table | TSV with unique `group` and all named finite environmental values |
| Geographic columns | Numerical `lon` and `lat` when `adjust_geography` is true |
| Permutation blocks | String `block` when `rda_blocked` is true. Define these from the sampling design before analysis. |

The genotype file may include extra individuals. The explicit sample table defines the cohort. Each ALT is a column. VCF haploid and diploid calls are supported, and group allele frequencies use the called ploidy denominator. Mixed-ploidy input should use VCF/BCF. The workflow handles the already-called panel, so upstream filtering and cohort eligibility are set before launch.

## Outputs

- `prepared/` includes aligned designs, explained variance, a locus eligibility catalogue and frozen method inputs.
- `lfmm_pc{3,5,10}/` gives coefficient, signed t statistic, nominal and GC-calibrated P, per-environment BH q, test status and calibration diagnostics.
- `rda_pc{3,5,10}/` gives global, marginal and axis tests, the permutation matrix, site scores, environmental arrows, locus loadings and exploratory candidates.
- `summary/` gives PC comparisons and descriptive overlap in clearly stated units.
- `run_fingerprint.json`, per-stage checksums and `progress.json` track provenance and resume status.

The default model has five joint environmental covariates, longitude and latitude, and LFMM2 K=5 held fixed. Free RDA permutations suit designs with exchangeable group residuals. For structured designs, evaluate blocked permutations before formal analysis. The primary PC specification, environmental panel and reporting thresholds should be set before comparing candidate sets.

This package provides generic workflow code, input templates and fictional examples. Place research data and formal outputs in local directories.
