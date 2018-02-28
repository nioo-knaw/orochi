__author__ = "Mattias de Hollander"
__copyright__ = "Copyright 2018, Mattias de Hollander"
__email__ = "m.dehollander@nioo.knaw.nl"
__license__ = "MIT"

from snakemake.utils import R

if os.path.isfile("config.json"):
    configfile: "config.json"

if True:
    include:
       "rules/data.input.rules"

if True:
    include:
        "rules/pre-processing/trimmomatic.rules"
if True:
    include:
        "rules/assembly/mapping.rules"
    include:
        "rules/assembly/stats.rules"

if config["assembler"] == "megahit":
    include:
        "rules/assembly/megahit.rules"


if config["assembler"] == "spades":
    include:
        "rules/assembly/spades.rules"

if True:
    include:
        "rules/binning/metabat.rules"

if True:
    include:
        "rules/report/report.rules"

output = []
output.append(rules.trimmomatic.output.fw_paired)
output.append(rules.report.output[0])
output.append(rules.metabat.output.depth)

rule final:
    input: expand(output,\
                  project=config["project"],\
                  sample=config["data"],\
                  treatment=config["treatment"],\
                  assembler=config["assembler"],\
                  kmers=config["assembly-klist"]\
                 )

