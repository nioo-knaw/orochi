__author__ = "Mattias de Hollander"
__copyright__ = "Copyright 2018, Mattias de Hollander"
__email__ = "m.dehollander@nioo.knaw.nl"
__license__ = "MIT"

from snakemake.utils import R

if os.path.isfile("config.json"):
    configfile: "config.json"

output = []

if True:
    include:
       "src/rules/data.input.smk"

if True:
    include:
#        "rules/pre-processing/trimmomatic.smk"
#    output.append(rules.trimmomatic.output.r2)
        "src/rules/pre-processing/bbduk.smk"
    output.append(rules.filter_contaminants.output[0])

# Read-based analysis
if False:
    include:
        "src/rules/read-based-analysis/fraggenescan.smk"
    output.append(rules.predict_genes.output[0])
    include:
        "src/rules/read-based-analysis/diamond.smk"
    output.append(rules.diamond_taxonomy_and_kegg.output.taxonomy)
    include:
        "src/rules/read-based-analysis/uproc.smk"
    output.append(rules.uproc.output.cog)
 
#if True:
#    include:
#        "src/rules/assembly/mapping.smk"
#    include:
#        "src/rules/assembly/stats.smk"

# Assembly
#if config["assembler"] == "megahit":
#    include:
#        "src/rules/assembly/megahit.smk"


#if config["assembler"] == "spades":
#    include:
#        "src/rules/assembly/spades.smk"

# Binning
#if True:
#    include:
#        "src/rules/binning/metabat.smk",
#    output.append(rules.metabat.output.depth)
#    include:
#        "src/rules/binning/mmgenome.smk",
#    output.append(rules.mmgenome_load_data.output[0])
     
# Reporting
#if True:
#    include:
#        "src/rules/report/report.smk"
#    output.append(rules.report.output[0])

rule final:
    input: expand(output,\
                  project=config["project"],\
                  sample=config["data"],\
                  treatment=config["treatment"],\
                  assembler=config["assembler"],\
                  kmers=config["assembly-klist"]\
                 )

