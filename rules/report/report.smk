rule create_rdata:
    input:
        quast="{project}/stats/quast.report.txt",
        flagstat="{project}/stats/flagstat.report.txt"
    output:
        rdata = "{project}/report/{project}.RData"
    run:
       R("""
       quast <- read.delim("{input.quast}")
       flagstat <- read.delim("{input.flagstat}")
       save.image(file="{output.rdata}")
       """)

rule report:
    input:
        rdata = "{project}/report/{project}.RData"
    output:
        "{project}/report/{project}.report.nb.html"
    params:
        prefix="{project}/report/{project}.report",
    conda: "../../envs/report.yaml"
    script:
        "report.Rmd"

