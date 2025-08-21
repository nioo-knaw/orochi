rule krona:
	input: rules.CAT.output.names
	output: 
		f"{outdir}/results/06_plots/{{sample_pool}}_krona.html"
	params:
		out_temp = temp(f"{outdir}/results/06_plots//{{sample_pool}}_contigs4krona_sep.txt")
	conda:
		"../envs/krona.yaml"
	threads:
		config['threads']
	shell:
		"""
		bash workflow/scripts/convert2krona.sh {input} > {params.out_temp}
		ktImportText {params.out_temp} -o {output}
		"""

