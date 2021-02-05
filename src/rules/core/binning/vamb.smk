rule sample_assembly:
    input: 
    output:
    conda: "../../../spades.yaml"
    shell:

rule concatenate:
    input:
    output:
    conda: "../../../vamb.yaml"
    shell:

rule read_mapper:
    input:
    output:
    conda: "../../../minimap2.yaml"
    shell:

rule vamb:
    input:
    output:
    conda: "../../../vamb.yaml"
    shell:

rule vamb_write_bins:
    input:
    output:
    conda: "../../../vamb.yaml"
    shell:
