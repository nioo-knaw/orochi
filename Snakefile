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
       "rules/data.input.smk"

if True:
    include:
        "rules/pre-processing/trimmomatic.smk"
    output.append(rules.trimmomatic.output.r2)

if True:
    include:
        "rules/assembly/mapping.smk"
    include:
        "rules/assembly/stats.smk"

# Assembly
if config["assembler"] == "megahit":
    include:
        "rules/assembly/megahit.smk"


if config["assembler"] == "spades":
    include:
        "rules/assembly/spades.smk"

# Binning
if True:
#    include:
#        "rules/binning/metabat.smk",
#    output.append(rules.metabat.output.depth)
    include:
        "rules/binning/mmgenome.smk",
    output.append(rules.mmgenome_load_data.output[0])
     
# Reporting
if True:
    include:
        "rules/report/report.smk"
    output.append(rules.report.output[0])

rule final:
    input: expand(output,\
                  project=config["project"],\
                  sample=config["data"],\
                  treatment=config["treatment"],\
                  assembler=config["assembler"],\
                  kmers=config["assembly-klist"]\
                 )

