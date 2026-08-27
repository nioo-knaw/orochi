# OrochiPipeline

Orochi is a Snakemake-based metagenomics pipeline designed for comprehensive analysis of metagenomic data, including preprocessing, assembly, binning, taxonomic and functional annotation, and Biosynthetic Gene Cluster (BGC) prediction.

The Orochi project is a collaboration between the Microbial Ecology group at the Netherlands Institute of Ecology (NIOO-KNAW) and the Bioinformatics Group at Wageningen University & Research (WUR).

## Overview

The Orochi pipeline integrates several state-of-the-art bioinformatics tools to process raw metagenomic reads and produce high-quality insights into microbial communities. It supports both single-sample assembly and co-assembly approaches.

### Key Features:
- **Preprocessing:** Quality control with `fastp` and host/contaminant removal with `bbmap`.
- **Assembly:** Supports `MEGAHIT` (for co-assembly) and `SPAdes` (for single-sample assembly).
- **Binning:** Multi-tool binning using `MetaBAT2` and `MaxBin2`, refined with `DAS Tool` and dereplicated with `dRep`.
- **Annotation:** 
    - Taxonomic profiling with `MetaPhlAn4` and `CAT`.
    - Functional annotation with `eggNOG-mapper`.
    - Eukaryotic gene prediction with `Augustus`.
- **BGC Prediction:** Identification of secondary metabolite clusters using `antiSMASH`.
- **MAG Linkage:** Linking MAGs using `markerMAG`.
- **Reporting:** Interactive HTML reports and visualizations.

## Requirements

- **Conda** or **Mamba** (recommended)
- **Snakemake** (>= 7.0.0)
- **Python** (>= 3.8)
- **R** (for reporting)

## Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/your-repo/orochi_nioo.git
   cd orochi_nioo
   ```

2. **Install Snakemake and Mamba:**
   It is highly recommended to use Mamba for environment management.
   ```bash
   conda install -c conda-forge mamba
   mamba create -c conda-forge -c bioconda -n snakemake snakemake
   ```

3. **Configure Databases:**
   The pipeline requires several external databases. Update the paths in `config/configfile.yaml`:
   - `phyloflash_db`
   - `emapper_database`
   - `antismash_db`
   - `CAT_database` & `CAT_taxonomy`
   - `checkm_db`

   For antiSMASH, you can set `download_antismash_db: true` to let the pipeline handle the installation.

## Usage

1. **Prepare the Sample Sheet:**
   Edit `config/samples.tsv`. It should be a tab-separated file with the following columns:
   - `sample`: Unique sample identifier.
   - `treatment1`, `treatment2`: Metadata columns.
   - `fq1`, `fq2`: Paths to raw paired-end FastQ files.
   - `sample_pool`: Grouping for co-assembly.

2. **Configure the Pipeline:**
   Edit `config/configfile.yaml` to set your output directory, assembly method, and tool parameters.

3. **Run the Pipeline:**
   Activate your Snakemake environment and run:
   ```bash
   snakemake --use-conda --cores <number_of_threads>
   ```

## Project Structure

```text
.
├── config/                 # Configuration files (configfile.yaml, samples.tsv)
├── resources/              # Reference genomes, contaminant sequences, and static data
├── workflow/
│   ├── envs/               # Conda environment definitions (.yaml)
│   ├── rules/              # Snakemake rule modules (.smk)
│   ├── scripts/            # Custom Python, R, and Bash scripts
│   └── Snakefile           # Main workflow entry point
├── README.md               # This file
└── LICENSE.md              # Project license
```

## Scripts and Workflow Modules

### Workflow Rules (`workflow/rules/`)
- `preprocessing.smk`: QC and filtering.
- `coassembly.smk` / `single_sample_assembly.smk`: Assembly logic.
- `binning.smk`: Binning, refinement, and dereplication.
- `prokaryote-gene_prediction_and_annotation.smk`: Taxonomic and functional profiling.
- `BGC_prediction.smk`: antiSMASH analysis.
- `plots.smk`: Report generation and visualization.

### Custom Scripts (`workflow/scripts/`)
- `render_report.R`: Renders the final HTML report.
- `augustify.py`: Wrapper for Augustus eukaryotic gene prediction.
- `summarize_antismash.py`: Aggregates antiSMASH results.
- `regenerate_antismash_html.py`: Utility for antiSMASH report formatting.

## Environment Variables

The pipeline primarily relies on the `configfile.yaml` for configuration. Ensure that any system-specific paths are correctly set there.

## Tests

TODO: Add instructions for running tests. A `test_rule.smk` exists in the workflow.

## License

TODO: Specify the license in `LICENSE.md`.
