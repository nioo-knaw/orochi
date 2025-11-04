rule krona:
	input: rules.CAT.output.names
	output: 
		f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_krona.html"
	params:
		out_temp = f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_contigs4krona_sep.txt"
	conda:
		"../envs/krona.yaml"
	threads:
		config['threads']
	shell:
		"""
		bash workflow/scripts/convert2krona.sh {input} > {params.out_temp}
		ktImportText {params.out_temp} -o {output}
		"""

rule bin_plots:
	input:
		checkm = f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample_pool}}_genomeinfo.tsv",
		bat = f"{outdir}/results/06_binning/BAT/{{sample_pool}}/{{sample_pool}}.bin2classification.names.txt"
		#bat = rules.BAT.output.bat_names
	output:
		tempfile = temp(f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_4scatterplot.tsv"),
		scatterplot = f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_bins_scatterplot.html",
		checkpoint = f"{outdir}/results/08_plots/.binplots_checkpoint/{{sample_pool}}.binPlots.done"
	params:
		tempfile = temp(f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_phyluminfo.tsv")
	threads:
		config['threads']
	resources:
		mem_mb=config['max_mem']
	log: f"{outdir}/logs/{{sample_pool}}_bins_scatterplot.log"
	conda:
		"../envs/python_clustering.yaml"
	shell:
		"""
		awk -v col=phylum 'NR==1{{for(i=1;i<=NF;i++){{if($i==col){{c=i;break}}}} print $c}} NR>1{{print $c}}' {input.bat} > {params.tempfile}
		paste -d',' {input.checkm} {params.tempfile} > {output.tempfile}
		python3 workflow/scripts/bin_scatterplot.py -i {output.tempfile} -o {output.scatterplot}
		touch {output.checkpoint}
		"""

SAMPLES_POOLS = glob_wildcards(f"{outdir}/results/06_binning/drep/checkm2_genomeinfo/{{sample_pool}}_genomeinfo.tsv").sample_pool
import glob

rule report:
    input:
        metaphlan_secondary = f"{outdir}/results/05_prokaryote_annotation/MetaPhlAn/merged_abundance_table.txt",
        binplots = expand(f"{outdir}/results/08_plots/{{sample_pool}}/{{sample_pool}}_bins_scatterplot.html", sample_pool=SAMPLES_POOLS)
        html_fastp = lambda wildcards: glob.glob(f"{outdir}/results/01_trimmed_reads/quality_reports/*.html")
    output:
        f"{outdir}/results/08_plots/Orochi_report.html"
    params:
        configfile= workflow.configfiles[0] if workflow.configfiles else "config/configfile.yaml",
        outdir_html = f"{outdir}/results/08_plots/rsc/"
    threads:
        config['threads']
    resources:
        mem_mb=config['max_mem']
    log: f"{outdir}/logs/report.log"
    conda:
        "../envs/html.yaml"
    shell:
        "mkdir {params.outdir_html}"
        "cp {input.html_fastp} {params.outdir_html}"
        "Rscript workflow/scripts/render_report.R {params.configfile} {input.metaphlan_secondary} {output}"
