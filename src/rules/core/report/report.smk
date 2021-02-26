"""
rule create_rdata:
    input:
        quast=expand("scratch/assembly/{assembler}/{treatment}/{kmers}/quast/report.txt", assembler=config["assembler"], treatment=config["treatment"], kmers=config["assembly-klist"]),
        flagstat="scratch/stats/flagstat.report.txt"
    output:
        rdata = "results/report/report.RData"
    run:
       R("""
       quast <- read.delim("{input.quast}")
       flagstat <- read.delim("{input.flagstat}")
       save.image(file="{output.rdata}")
       """)

rule report:
    input:
        rdata = "results/report/report.RData"
    output:
        "results/report/report.nb.html"
    params:
        prefix="results/report/report",
    conda: "../../../envs/report.yaml"
    script:
        "report.Rmd"
"""
