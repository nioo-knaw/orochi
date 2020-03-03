__author__ = "Mattias de Hollander"
__copyright__ = "Copyright 2018, Mattias de Hollander"
__email__ = "m.dehollander@nioo.knaw.nl"
__license__ = "MIT"

from snakemake.utils import R, listfiles

if os.path.isfile("config.json"):
    configfile: "config.json"


# List all modules/rules
modules = list(listfiles('src/rules/{module}/{mod}.smk'))
for mod in modules:
    print(mod)

if True:
    include:
       "src/rules/data.input.smk"

if True:
    include:
#        "rules/pre-processing/trimmomatic.smk"
        "src/rules/pre-processing/bbduk.smk"

# Read-based analysis
if False:
    include:
        "src/rules/read-based-analysis/fraggenescan.smk"
    include:
        "src/rules/read-based-analysis/diamond.smk"
    include:
        "src/rules/read-based-analysis/uproc.smk"
 
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
#    include:
#        "src/rules/binning/mmgenome.smk",
     
# Reporting
#if True:
#    include:
#        "src/rules/report/report.smk"

# Dynamically add all output files
output = []
for name,rule in rules.__dict__.items():
    for file in rule.output:
        output.append(file)

rule final:
    input: expand(output,\
                  project=config["project"],\
                  sample=config["data"],\
                  treatment=config["treatment"],\
                  assembler=config["assembler"],\
                  kmers=config["assembly-klist"]\
                 )

