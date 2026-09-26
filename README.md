# PMM2-CDG Multiomic Integration and Robustness Workflow

This repository contains the computational workflow developed for the study
of multiomic integration in PMM2-CDG.

The workflow integrates phenotypic, transcriptomic, miRNA, metabolomic,
genetic and clinical-severity information to investigate molecular patterns
associated with disease severity. In addition to the main multiomic analysis,
the repository includes complementary sensitivity and robustness analyses
designed to assess the influence of data-layer selection, individual patients,
individual network features and stochastic network embedding.

The workflow is intended for exploratory and hypothesis-generating analyses
in the context of small rare-disease cohorts.

---

## Workflow overview

The analysis is organized into the following stages:

1. Environment setup
2. STRING interaction-network retrieval and preprocessing
3. Dataset parsing and preprocessing
4. Multi-Factor Analysis (MFA)
5. Leave-one-layer-out sensitivity analysis
6. Leave-one-sample-out robustness analysis
7. Posterior network and functional analyses
8. Node2vec/HDBSCAN stability analysis
9. Workflow monitoring and recovery

---

## Software environment

The analyses were performed using Python and R.

### Python environment

- Python 3.14.3
- NumPy 2.4.0
- pandas 2.3.3
- SciPy 1.17.1
- scikit-learn 1.8.0
  - Elastic-net modelling
  - `sklearn.cluster.HDBSCAN`
- scikit-bio 0.7.0
  - PERMANOVA
  - PERMDISP
- statsmodels 0.14.6
- PecanPy 2.0.9
  - node2vec network embedding

Workflow-specific software:

- PETS (`py_pets`) 1.3.6
- `py_exp_calc` 1.2.1
- AutoFlow Next 1.1.7
- NetAnalyzer: 1.1.0

### R environment

- R 4.5.1
- Bioconductor 3.22

The complete R `sessionInfo()` used to document the computational environment
is provided in the Supplementary Material associated with the manuscript.

---

## Reproducibility settings

Random seeds are fixed for stochastic procedures whenever supported by the
underlying software.

The principal workflow-level seed is:

```text
seed = 123
```

Additional stochastic settings and seeds used for individual analyses are
defined in the workflow templates and scripts included in this repository.

Important analysis parameters include:

```text
MFA / workflow seed: 123
MFA-associated feature cutoff: |r| > 0.70
STRING confidence cutoff: 700
miRNA interaction database: miRecords
```

For the network embedding analysis, the node2vec and HDBSCAN parameters used
in the study are documented in the corresponding workflow templates.

The node2vec stability analysis repeats the complete embedding and clustering
procedure using multiple random seeds while keeping all other hyperparameters
fixed.

---

## Data availability and directory structure

Patient-level multiomic and clinical datasets are not publicly distributed
because they contain sensitive information from individuals affected by an
ultra-rare disease and are subject to ethical and data-protection restrictions.

A private repository containing the patient datasets is associated with this
project. Access requires explicit authorization.

The workflow expects the private dataset repository to be available under:

```text
datasets/
```

The analysis script currently defines:

```text
datasets/datasets/
```

as the source-data directory, while parsed datasets are written to:

```text
parsed_datasets/
```

The expected project structure is therefore approximately:

```text
project/
├── workflow.sh
├── aux_scripts/
├── templates/
├── datasets/
│   └── datasets/
│       └── sample_cohort/
└── parsed_datasets/
```

The `sample_cohort` directory contains the data layers required by the
multiomic workflow, including:

```text
phenotypes.txt
all_genes_var_filtered.txt
label_genes.txt
miRNA.txt
metabolomics.txt
severity.txt
severity_scales.txt
variants.txt
```

The exact patient-level files are not distributed publicly.

---

# Workflow execution

The workflow is controlled through a single shell script:

```bash
./workflow.sh <mode> [options]
```

The execution modes currently implemented are described below.

---

## 1. Environment setup

```bash
./workflow.sh set_env
```

This mode initializes the local Python virtual environment used by the
workflow.

It should be run before the analysis when a local environment has not already
been configured.

---

## 2. Retrieve and preprocess the STRING interaction network

```bash
./workflow.sh get_string
```

This step downloads the human STRING v12.0 protein-protein interaction network
and corresponding identifier mappings.

Interactions are filtered using:

```text
STRING combined score >= 700
```

The resulting network is converted to Ensembl- and HGNC-compatible
representations for downstream network analysis.

---

## 3. Dataset parsing and preprocessing

```bash
./workflow.sh parse
```

This stage prepares the sample-cohort datasets used by the subsequent
analyses.

The parsing procedure includes:

- phenotype/HPO annotation;
- gene-identifier annotation;
- stability-based mRNA feature selection;
- miRNA annotation;
- preparation of metabolomic data;
- preparation of severity variables;
- preparation of clinical severity scales;
- preparation of genetic-variant information.

The resulting processed datasets are written to:

```text
parsed_datasets/sample_cohort/
```

The mRNA stability-selection procedure is executed using:

```bash
stable_select
```

with the parameters documented in the manuscript and associated software.

---

## 4. Multi-Factor Analysis

```bash
./workflow.sh ma [AutoFlow options]
```

This step runs the primary multiomic integration using AutoFlow.

The MFA integrates the active molecular and phenotypic data layers and
projects supplementary clinical variables without allowing them to contribute
to construction of the latent dimensions.

The workflow-level random seed is:

```text
123
```

Execution results are stored under the configured MFA execution directory.

---

## 5. Leave-one-layer-out sensitivity analysis

```bash
./workflow.sh select_layers [AutoFlow options]
```

This mode evaluates the contribution of each active data layer to the MFA
integration.

Each active layer is independently removed and the MFA is recomputed using
the remaining layers.

The resulting configurations are compared with the complete integration using:

- RV coefficient;
- RV2 coefficient;
- absolute change in severity-associated PERMANOVA effect size,
  `|Delta R2|`.

This analysis was used to assess the contribution of:

- gene expression;
- miRNA expression;
- phenotypic variables;
- metabolomics.

---

## 6. Leave-one-sample-out robustness analysis

```bash
./workflow.sh check_outlier_robustness [AutoFlow options]
```

This mode assesses whether the feature-selection and MFA results are
disproportionately influenced by individual patients.

Each patient is removed once and the relevant analysis steps are repeated.

Two levels of robustness are evaluated.

### Feature-selection robustness

The stability-selected gene set obtained from each leave-one-sample-out run
is compared with the full-cohort reference set using:

- Jaccard similarity;
- reference-set retention.

### MFA-associated feature robustness

The MFA is recomputed after removing each patient.

MFA dimensions are matched to the corresponding full-cohort dimensions, and
genes satisfying:

```text
|r| > 0.70
```

are compared with the corresponding full-cohort PC1- and PC2-associated gene
sets.

Jaccard similarity and reference-set retention are used to quantify
agreement.

---

## 7. Posterior network and functional analysis

```bash
./workflow.sh pa [AutoFlow options]
```

This stage performs the downstream analyses based on the MFA results.

The posterior workflow includes:

- selection of features strongly associated with MFA dimensions;
- mapping of selected genes onto the protein-protein interaction network;
- node2vec network embedding;
- HDBSCAN module detection;
- module-level PERMANOVA;
- PERMDISP assessment;
- leave-one-feature-out module sensitivity analysis;
- functional enrichment;
- miRNA-module integration.

The cutoff used to select genes associated with MFA dimensions is:

```text
|r| > 0.70
```

The configured miRNA interaction resource is:

```text
miRecords
```

### Exact PERMANOVA testing

For the five High- and five Low-severity individuals, module-level PERMANOVA
uses exhaustive enumeration of all:

```text
252
```

possible balanced 5-versus-5 group assignments.

For each module, the analysis reports:

- PERMANOVA pseudo-F;
- PERMANOVA R2;
- exact raw P value;
- Benjamini-Hochberg FDR-adjusted P value.

### PERMDISP

PERMDISP is applied to the same distance matrices to assess whether significant
PERMANOVA results may be influenced by differences in within-group
multivariate dispersion.

The same exhaustive set of 252 balanced assignments is used.

### Leave-one-feature-out sensitivity analysis

For PC1-derived modules showing significant PERMANOVA associations, each
network feature is individually removed while module membership is kept fixed.

PERMANOVA is then recomputed on the remaining features.

Feature influence is quantified as:

```text
Delta R2 = R2_full - R2_LOFO
```

This analysis evaluates whether a module-level severity association is
disproportionately dependent on an individual feature.

---

## 8. Node2vec/HDBSCAN stability analysis

```bash
./workflow.sh check_embedding_stability [AutoFlow options]
```

This mode assesses the stochastic robustness of network-module identification.

The complete node2vec embedding and HDBSCAN clustering procedure is repeated
using distinct random seeds for the PC1- and PC2-associated gene sets while
keeping embedding and clustering hyperparameters fixed.

The analysis comprises 100 repeated embedding/clustering realizations per
dimension.

Each resulting partition is compared with the reference partition using:

- Adjusted Rand index (ARI);
- Normalized mutual information (NMI);
- co-assignment Jaccard;
- best-match Jaccard.

ARI and NMI quantify global similarity between partitions.

Co-assignment Jaccard evaluates preservation of pairwise feature
co-clustering.

Best-match Jaccard quantifies correspondence between individual reference
modules and their closest matching modules across repeated runs.

---

## 9. Workflow status

```bash
./workflow.sh check
```

This command reports the status of the main MFA workflow through the AutoFlow
execution log.

---

## 10. Workflow recovery

```bash
./workflow.sh recover
```

This command can be used to inspect and recover interrupted AutoFlow analyses.

---

# Analysis parameters

The principal parameters exposed by the top-level workflow include:

```text
STRING confidence threshold = 700
Global workflow seed = 123
MFA-associated feature cutoff = |r| > 0.70
miRNA interaction resource = miRecords
```

Additional method-specific parameters are defined in the scripts and AutoFlow
templates distributed with this repository.

These include parameters for:

- stability selection;
- MFA;
- node2vec;
- HDBSCAN;
- PERMANOVA/PERMDISP;
- leave-one-layer-out analysis;
- leave-one-sample-out analysis;
- leave-one-feature-out analysis;
- repeated embedding stability assessment.

---

# Reproducing the analysis

A complete analysis is typically executed in the following order:

```bash
./workflow.sh set_env
./workflow.sh get_string
./workflow.sh parse
./workflow.sh ma
./workflow.sh select_layers
./workflow.sh check_outlier_robustness
./workflow.sh pa
./workflow.sh check_embedding_stability
```

Execution status can be inspected at any time with:

```bash
./workflow.sh check
```

and interrupted workflows can be inspected/recovered with:

```bash
./workflow.sh recover
```

Because patient-level data are restricted, full reproduction of the biological
results requires authorized access to the private dataset repository.

The publicly available repository nevertheless provides the analysis code,
workflow definitions, computational parameters and software environment
information required to reproduce the computational procedure when compatible
input data are available.

---

# Notes

- Patient-level datasets must not be committed to the public repository.
- Parsed datasets and intermediate analysis outputs should remain untracked.
- The workflow is intended for exploratory and hypothesis-generating analyses.
- Robustness analyses provide internal sensitivity assessments and do not
  substitute validation in an independent patient cohort.
- Exact software versions should be preserved when reproducing the reported
  analysis.