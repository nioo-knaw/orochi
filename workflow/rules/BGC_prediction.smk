""" The rules related to Biosynthetic Gene Cluster (BGC) prediction and related analyses"""

rule antismash:
    input:
        "path/to/contigs"

    output:
        "path/to/output"

    shell:
        "antismash_meta -options"

rule bigscape:
    input:
        "path/to/antismash_output"

    output:
        "path/to/output"

    shell:
        "bigscape -options"

rule itol_bgc:
    input:
        "path/to/input"
    output:
        "path/to/output"
