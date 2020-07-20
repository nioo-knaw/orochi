**Note: This is the development version of Orochi, use at own risk.**

## Orochi 1.0 (alpha)
This repo is a rewrite of the initial Orochi pipeline for shotgun metagenomics and still in development. It does not perform a full analysis yet.

### What changed compared to Orochi 0.1?
  - Rules are now organized into modules. 
    - There is a set of core modules that do the QC, assembly and annotation.
    - Additional modules can be added as extension, for example the detection of secondary metabolite biosynthesis gene clusters with antismash
    - Modules are automatically detected in the main Snakefile
  - The QC step is now performed with bbduk in stead of trimmomatic to 1) trim poly-G tails from Novaseq sequences 2) remove low complexity reads 3) remove phix reads

### Getting started

1. Logon to the place where you will analysis your data, e.g. server
2. Create a local copy of the pipeline in a project folder
```
git clone https://gitlab.bioinf.nioo.knaw.nl/pipelines/orochi -b modules
``` 
3. Enter your (NIOO) username and password

This will create a local copy of the pipeline in a folder called `orochi`

4. Enter the pipeline folder with `cd orochi`

5. The configuration of the pipeline needs to be set in a file called `config.json`. There is an example file which you can copy and adjust
```
cp config.sample.json config.json
```

Take a look at the `config.json` file

6. Adjust settings like:
   - project name
   - host genome
   - data files
   - treatment/assembly groups

7. Check if the config file is correct and which steps will be run
```
snakemake -n
```

8. Look at a diagram of all the steps in the workflow

```
snakemake --rulegraph | dot -Tpng > workflow.png
display workflow.png
```

9. Run the pipeline. `-j` specifies the number of threads. Conda is the package manager. Optionally do this in a tmux session.
```
snakemake -j 8 --use-conda
```