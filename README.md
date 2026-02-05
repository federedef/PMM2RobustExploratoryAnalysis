# PMM2 Robust Exploratory Analysis Workflow

This repository contains a computational workflow developed for the article on PMM2-related disorders.
The goal of this workflow is to perform a robust exploratory analysis on multi-omics and clinical data
from patients with PMM2 disease, with the aim of identifying potential biomarkers that discriminate
between patients with high severity and low severity phenotypes.

The workflow integrates phenotypic, genetic, transcriptomic, miRNA, and metabolomic data and applies
a multifactorial analysis framework followed by posterior analyses to explore severity-associated
patterns in the data.

--------------------------------------------------
OVERVIEW
--------------------------------------------------

The workflow is organized into four main stages:

1. Environment setup
2. Dataset parsing and preprocessing
3. Multifactor (robust exploratory) analysis
4. Posterior analysis and result inspection

Parsed datasets and intermediate results are generated locally and are not intended to be versioned.

--------------------------------------------------
REQUIREMENTS
--------------------------------------------------

- Bash-compatible shell
- Python (with venv support)
- R (with required annotation scripts available in the environment)
- External tools available in the environment:
  - semtools
  - stable_select
  - AutoFlow

Environment-specific initialization scripts are assumed to be available and correctly configured.

--------------------------------------------------
DATA STRUCTURE
--------------------------------------------------

The input datasets required to run this workflow are not publicly distributed due to data privacy constraints.

A separate **private repository** containing the patient datasets is associated with this project and is
integrated into the main repository as a **git submodule**. Users must have explicit access permissions
to this private repository in order to run the workflow.

Once access is granted, the private dataset repository is expected to be available locally under a
folder named `datasets`.

Parsed and intermediate datasets generated during execution are created locally in a separate input
folder and are not intended to be versioned or shared.

--------------------------------------------------
WORKFLOW EXECUTION
--------------------------------------------------

The workflow is controlled via a single script with different execution modes.

1) Environment setup

This step creates a Python virtual environment and installs required dependencies.

Command:
./workflow.sh set_env

Note:
This step must be executed once before running the parsing stage.

--------------------------------------------------

2) Dataset parsing and preprocessing

This step parses and preprocesses all datasets, generating normalized and annotated versions
of phenotypes, genes, miRNAs, and metabolomic data for both full and sample cohorts.

Command:
./workflow.sh parse

Requirements:
- The virtual environment must exist
- The `datasets` folder must be present and correctly populated

--------------------------------------------------

3) Multifactor analysis (robust exploratory analysis)

This step performs the core multifactorial exploratory analysis aimed at identifying
patterns associated with disease severity.

Command:
./workflow.sh ma [options]

This stage uses AutoFlow to execute the multifactor analysis workflow.

--------------------------------------------------

4) Posterior analysis

This step performs posterior analyses based on the results of the multifactor analysis,
allowing further interpretation and refinement of severity-associated signals.

Command:
./workflow.sh pa [options]

--------------------------------------------------

5) Workflow monitoring and recovery

Check the status of running or completed analyses:
./workflow.sh check

Recover or resume interrupted analyses:
./workflow.sh recover

--------------------------------------------------
NOTES
--------------------------------------------------

- Parsed datasets are generated locally and should not be committed to the repository.
- The workflow is designed for exploratory and hypothesis-generating analyses.
- Results should be interpreted in the context of PMM2 disease biology and validated independently.

--------------------------------------------------