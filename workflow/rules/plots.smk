from workflow.scripts.augustify import outdir

rule krona:
    input: rules.CAT.output.names
    output:
        f"{outdir}/results/09_plots/rsc/{{sample_pool}}/{{sample_pool}}_krona.html"
    params:
        out_temp = f"{outdir}/results/09_plots/rsc/{{sample_pool}}/{{sample_pool}}_contigs4krona_sep.txt",
        temp_dir = f"{outdir}/results/09_plots/rsc/{{sample_pool}}"
    conda:
        "../envs/krona.yaml"
    threads:
        config['threads']
    log:
        f"{outdir}/logs/krona/krona_{{sample_pool}}.log"
    shell:
        """
        mkdir -p {params.temp_dir}
        bash workflow/scripts/convert2krona.sh {input} > {params.out_temp} 2>> {log}
        ktImportText {params.out_temp} -o {output} >> {log} 2>&1
        echo "Krona plot successfully generated for {wildcards.sample_pool}" >> {log}
        """

# SAMPLES_POOLS = glob_wildcards(f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample_pool}}_genomeinfo.tsv").sample_pool
import glob
import json

rule report:
    input:
        metaphlan_secondary = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.txt",
        antismash_bac = expand(f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial/index.html",
            sample_pool=ASSEMBLY_UNITS
        ),
        antismash_fun = expand(f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal/index.html",
            sample_pool=ASSEMBLY_UNITS
        ),
        html_fastp = expand(f"{outdir}/results/01_trimmed_reads/quality_reports/{{sample}}.html",
            sample=samples["sample"]
        ),
        markermag_done = expand(f"{outdir}/results/07_maglinkage/{{sample_pool}}/markermag/markermag.done",
            sample_pool=ASSEMBLY_UNITS)
    output:
        f"{outdir}/results/09_plots/Orochi_report.html"
    params:
        configfile= workflow.configfiles[0] if workflow.configfiles else "config/configfile.yaml",
        outdir_html = f"{outdir}/results/09_plots/rsc/",
        antismash_bac = json.dumps([
            {
                "src": f"{outdir}/results/08_BGC/antismash/{sp}/bacterial",
                "dst": f"{outdir}/results/09_plots/rsc/{sp}/antismash_bac/"
            }
            for sp in ASSEMBLY_UNITS
        ]),
        antismash_fun = json.dumps([
            {
                "src": f"{outdir}/results/08_BGC/antismash/{sp}/fungal",
                "dst": f"{outdir}/results/09_plots/rsc/{sp}/antismash_fun/"
            }
            for sp in ASSEMBLY_UNITS
        ]),
        rep_antismash_bac = expand(f"{outdir}/results/09_plots/rsc/{{sample_pool}}/antismash_bac/",
            sample_pool=ASSEMBLY_UNITS
        ),
        rep_antismash_fun = expand(f"{outdir}/results/09_plots/rsc/{{sample_pool}}/antismash_fun/",
            sample_pool=ASSEMBLY_UNITS
        )
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    log: f"{outdir}/logs/report.log"
    conda:
        "../envs/html.yaml"
    shell:
        r"""
        set -euo pipefail
        
        echo "Starting Orochi report generation..." > {log}
        echo "Timestamp: $(date)" >> {log}
        
        mkdir -p {params.outdir_html}
        echo "Copying fastp HTML reports..." >> {log}
        cp {input.html_fastp} {params.outdir_html} 2>> {log}
        
        mkdir -p {params.rep_antismash_bac}
        mkdir -p {params.rep_antismash_fun}
        
        echo "Copying bacterial antiSMASH results..." >> {log}
        echo '{params.antismash_bac}' | jq -c '.[]' | while read pair; do
            src=$(echo "$pair" | jq -r '.src')
            dst=$(echo "$pair" | jq -r '.dst')
            echo "  Processing: $src -> $dst" >> {log}
            mkdir -p "$dst"
            find "$src" -type f \
                \( -name "*.html" -o -name "*.js" -o -name "*.css" -o -name "*.svg" -o -name "*.png" \) \
                -print0 | while IFS= read -r -d '' file; do
                    rel=${{file#"$src"/}}
                    mkdir -p "$dst/$(dirname "$rel")"
                    cp "$file" "$dst/$rel"
                done
        done 2>> {log}
        
        echo "Copying fungal antiSMASH results..." >> {log}
        echo '{params.antismash_fun}' | jq -c '.[]' | while read pair; do
            src=$(echo "$pair" | jq -r '.src')
            dst=$(echo "$pair" | jq -r '.dst')
            echo "  Processing: $src -> $dst" >> {log}
            mkdir -p "$dst"
            find "$src" -type f \
                \( -name "*.html" -o -name "*.js" -o -name "*.css" -o -name "*.svg" -o -name "*.png" \) \
                -print0 | while IFS= read -r -d '' file; do
                    rel=${{file#"$src"/}}
                    mkdir -p "$dst/$(dirname "$rel")"
                    cp "$file" "$dst/$rel"
                done
        done 2>> {log}
        
        echo "Rendering final HTML report with R..." >> {log}
        
        Rscript workflow/scripts/render_report.R {params.configfile} {input.metaphlan_secondary} {output} >> {log} 2>&1
        
        echo "Report generation completed successfully at $(date)" >> {log}
        """
