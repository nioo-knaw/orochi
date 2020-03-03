rule create_rdata:
    input:
        quast="scratch/stats/quast.report.txt",
        flagstat="scratch/stats/flagstat.report.txt"
    output:
        rdata = "scratch/report/orochi.RData"
    run:
       R("""
       quast <- read.delim("{input.quast}")
       flagstat <- read.delim("{input.flagstat}")
       save.image(file="{output.rdata}")
       """)

rule report:
    input:
        rdata = "scratch/report/orochi.RData",
        mmgenome=expand("scratch/binning/mmgenome/{assembler}/{treatment}/{kmers}/orochi.RData", assembler=config["assembler"], treatment=config["treatment"], kmers=config["assembly-klist"])
    output:
        "scratch/report/scratch.report.nb.html"
    params:
        prefix="scratch/report/scratch.report",
    conda: "../../envs/report.yaml"
    script:
        "report.Rmd"

