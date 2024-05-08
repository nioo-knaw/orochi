# OrochiPipeline


##DRY RUN
snakemake -npc1 --use-conda --conda-frontend conda results/02_filtered_reads/B55_filt_1.fastq.gz

##RUN
snakemake -pc16 --use-conda --conda-frontend conda results/02_filtered_reads/B55_filt_1.fastq.gz
