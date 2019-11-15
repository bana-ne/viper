#!/usr/bin/env python

def optitype_targets(wildcards):
    """Generates the targets for this module"""
    ls = []
    for s in config["ordered_sample_list"]:
        ls.append("analysis/optitype/%s/%s_result.tsv" % (s,s))
        ls.append("analysis/optitype/%s/%s_coverage_plot.pdf" % (s,s))
    return ls

rule optitype_all:
    input:
        optitype_targets
        
rule optitype_get_chr6_reads:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam",
        bai="analysis/STAR/{sample}/{sample}.sorted.bam.bai"
    output:
        "analysis/STAR/{sample}/{sample}.sorted.chr6.bam"
    message: "OPTITYPE: obtaining chr6 reads"
    shell:
        "samtools view -b -h {input.bam} chr6 > {output}"

rule optitype_convert_bam2fq:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.chr6.bam"
    output:
        #USING _mates from align.snakefile
        mates = expand( "analysis/STAR/{{sample}}/{{sample}}.chr6.{mate}.fq", mate=_mates)
    message: "OPTITYPE: converting chr6 bam to fastq"
    run:
        #NOTE: need to run things differently if SE vs PE
        if len(_mates) == 1: #SE--redirect output
            shell("samtools bam2fq {input} > {output.mates}")
        else: #PE --use the -1 and -2 params
            out_files = " ".join(["-%s %s" % (i+1,m) for (i,m) in enumerate(output.mates)])
            shell("samtools bam2fq {out_files} {input}")

rule optitype:
    input:
        mates = expand( "analysis/STAR/{{sample}}/{{sample}}.chr6.{mate}.fq", mate=_mates)
    output:
        results="analysis/optitype/{sample}/{sample}_result.tsv",
        coverage="analysis/optitype/{sample}/{sample}_coverage_plot.pdf"
    message: "OPTITYPE: running optitype"
    params:
        pypath="PYTHONPATH=%s" % config["python2_pythonpath"],
        outpath="analysis/optitype/{sample}",
        name="{sample}",
        python2=config['python2'],
        optitype_path=config['optitype_path'],
    shell:
        "{params.pypath} {params.python2} {params.optitype_path}/OptiTypePipeline.py -i {input} -r -o {params.outpath} -c viper/static/optitype/config.ini -p {params.name}"
