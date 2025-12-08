rule krona:
	input: rules.CAT.output.names
	output: 
		f"{outdir}/results/09_plots/{{sample_pool}}/{{sample_pool}}_krona.html"
	params:
		out_temp = f"{outdir}/results/09_plots/{{sample_pool}}/{{sample_pool}}_contigs4krona_sep.txt"
	conda:
		"../envs/krona.yaml"
	threads:
		config['threads']
	shell:
		"""
		bash workflow/scripts/convert2krona.sh {input} > {params.out_temp}
		ktImportText {params.out_temp} -o {output}
		"""

SAMPLES_POOLS = glob_wildcards(f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample_pool}}_genomeinfo.tsv").sample_pool
import glob
import json

rule report:
    input:
        metaphlan_secondary = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.txt",
        antismash_bac = expand(f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/bacterial", sample_pool=SAMPLES_POOLS),
        antismash_fun = expand(f"{outdir}/results/08_BGC/antismash/{{sample_pool}}/fungal", sample_pool=SAMPLES_POOLS),
        html_fastp = lambda wildcards: glob.glob(f"{outdir}/results/01_trimmed_reads/quality_reports/*.html")
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
            for sp in SAMPLES_POOLS
        ]),
        antismash_fun = json.dumps([
            {
                "src": f"{outdir}/results/08_BGC/antismash/{sp}/fungal",
                "dst": f"{outdir}/results/09_plots/rsc/{sp}/antismash_fun/"
            }
            for sp in SAMPLES_POOLS
        ]),
        rep_antismash_bac = expand(f"{outdir}/results/09_plots/rsc/{{sample_pool}}/antismash_bac/", sample_pool=SAMPLES_POOLS),
        rep_antismash_fun = expand(f"{outdir}/results/09_plots/rsc/{{sample_pool}}/antismash_fun/", sample_pool=SAMPLES_POOLS)
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    log: f"{outdir}/logs/report.log"
    conda:
        "../envs/html.yaml"
    shell:
        """
        mkdir -p {params.outdir_html}
        cp {input.html_fastp} {params.outdir_html}
        mkdir -p {params.rep_antismash_bac}
        mkdir -p {params.rep_antismash_fun}
        echo '{params.antismash_bac}' | jq -c '.[]' | while read pair; do
            src=$(echo "$pair" | jq -r '.src')
            dst=$(echo "$pair" | jq -r '.dst')
            mkdir -p "$dst"
            find "$src" -type f \
                \( -name "*.html" -o -name "*.js" -o -name "*.css" -o -name "*.svg" -o -name "*.png" \) \
                -print0 | while IFS= read -r -d '' file; do
                    rel=${{file#"$src"/}}
                    mkdir -p "$dst/$(dirname "$rel")"
                    cp "$file" "$dst/$rel"
                done
        done
        echo '{params.antismash_fun}' | jq -c '.[]' | while read pair; do
            src=$(echo "$pair" | jq -r '.src')
            dst=$(echo "$pair" | jq -r '.dst')
            mkdir -p "$dst"
            find "$src" -type f \
                \( -name "*.html" -o -name "*.js" -o -name "*.css" -o -name "*.svg" -o -name "*.png" \) \
                -print0 | while IFS= read -r -d '' file; do
                    rel=${{file#"$src"/}}
                    mkdir -p "$dst/$(dirname "$rel")"
                    cp "$file" "$dst/$rel"
                done
        done
        Rscript workflow/scripts/render_report.R {params.configfile} {input.metaphlan_secondary} {output}
        """
