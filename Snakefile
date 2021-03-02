__author__ = "Mattias de Hollander"
__copyright__ = "Copyright 2020, Mattias de Hollander"
__email__ = "m.dehollander@nioo.knaw.nl"
__license__ = "GNU GPLv3"

from snakemake.utils import R, listfiles

if os.path.isfile("config.json"):
    configfile: "config.json"

# Dynamically load all modules/rules
smks = list(listfiles('src/rules/core/{section}/{part}.smk'))

for smk,rule in smks:
    include: smk

# Load extensions
"""
for ext, rule in list(listfiles('src/rules/extensions/{section}/{part}.smk')):
    section, name = rule
    if section in ["antismash"]:
        include: ext
"""

#Load report
"""
for rpt, rule in list(listfiles('src/rules/report/{part}.smk)):
    name = rule
    include: rpt
"""

# Dynamically add all output files
output = []
for name,rule in rules.__dict__.items():
    for file in rule.output:
        output.append(file)

rule final:
    input: expand(output,\
                  sample=config["data"],\
                  treatment=config["treatment"],\
                  assembler=config["assembler"],\
                  binner=config["binner"],\
                  kmers=config["assembly-klist"]\
                 )

