# N4M Results Summry

**Nextflow pipeline for summarizing results obtained from other data processing workflows in Nextflow4Metabolomics Suite**.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->
**N4M Results Summry** is a workflow that summarize LC-HRMS metabolomics data processing results obtained from other workflows in the Nextflow4Metabolomics Suite.

- The pipeline containers two processes: 1) merging peaks from two peak tables, the two peak tables can be obtained by runnin the same workflow using two OS; 2) quantifiying peaks before and after merging.

- The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

- The workflow support MacOS and Linux operating systems.

## Installation

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility.

## Quick Start

1. Download the pipeline repo and dirct to the folder:
	```bash
	git clone https://github.com/Nextflow4Metabolomics/n4m_results_summary.git && cd n4m_results_summary
	```

2. Run the pipeline with example data:
    ```bash
    nextflow main.nf -profile functional_test > logs/execution.log
    ```

## Process Your Own Data

1. Put the data files in the folder `data/`, and named the two files `table_1` and `table_2`, note that `table_2` should have less peaks than `table_1`.
2. Run the pipeline (use "docker" as the profile when running locally, and "singularity" as the profile when running with a high-performance computing system):
    ```bash
    nextflow main.nf -profile docker > logs/execution.log
    ```

## Configuration

- Configuration for running with Docker are set in the file `conf/base.config`.
- Configuration for running with High-Performance Computing and Singularity are set in the file `conf/HiPerGator.config`.
- Parameters for MS-DIAL and MS-FLO are set in their specific configuration files.

## Credits

Dr. Dominick Lemas (Xinsong Du's Ph.D. advisor) and Xinsong Du play an important role on conceptulization.
The `nextflow4ms-dial` was mainly developed by Xinsong Du. 

We thank the following people for their extensive assistance in the development
of this pipeline:
