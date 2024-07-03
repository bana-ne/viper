#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab

#-------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-------------------------------
# =============================================================================
# Copyright: Copyright 2024, viper_adaption
# License: GPL
# Version: v0.1
# Editor: Vanessa Mandel
# Email: schmollv@mytum.de
# Created: 27.05.2024 (DD.MM.YYYY)
# Last Modified: 
# =============================================================================
""" Description:
    Rules specific to STAR and STAR-Fusion (identifies potential fusion genes discoered during the alignment) from the Viper pipeline
    Cornwell M, et al. VIPER: Visualization Pipeline for RNA-seq, a Snakemake workflow for efficient and complete RNA-seq analysis. BMC Bioinformatics. 2018 PMID: 29649993.

    The script was adapted to optionally use previously aligned files from other pipelines.
"""
# =============================================================================

# config["samples"] contains all sample names defined in the config.yaml file under samples
def getFastq(wildcards):
    return config["samples"][wildcards.sample]
# get paired-end file names for fusion:
def getLeftFastq(wildcards):
    return config["samples"][wildcards.sample][0]
def getRightFastq(wildcards):
    return config["samples"][wildcards.sample][1]    


strand_command=""
rRNA_strand_command=""

# set strand specific filtering options for STAR alignment
if( config["stranded"] ):
    strand_command="--outFilterIntronMotifs RemoveNoncanonical"
    rRNA_strand_command="--outFilterIntronMotifs RemoveNoncanonical"
else:
    rRNA_strand_command="--outSAMstrandField intronMotif"

run_fusion= True if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else False # turn run_fusion to true if paired-end reads are available
gz_command="--readFilesCommand zcat" if config["samples"][config["ordered_sample_list"][0]][0][-3:] == '.gz' else "" # check if input files need to be decompressed

if( run_fusion ):
    strand_command = """ --outFilterIntronMotifs RemoveNoncanonicalUnannotated --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 200000 --alignIntronMax 200000 """
_mates = ['mate1', 'mate2'] if len(config["samples"][config["ordered_sample_list"][0]]) == 2 else ['mate1']
_keepPairs = "KeepPairs" if len(_mates) == 2 else ""

# definition of the STAR alignment process
rule run_STAR:
    """ STAR alignment command using parameters that later allow STAR-Fusion to run on the output 
        (i.e. using the chimeric read detection & reporting in Chimeric.junction.out files)
    """
    input:
        getFastq
    output:
        bam=protected("analysis/STAR/{sample}/{sample}.sorted.bam"),
        counts="analysis/STAR/{sample}/{sample}.counts.tab",
        log_file="analysis/STAR/{sample}/{sample}.Log.final.out",
        #COOL hack: {{sample}} is LEFT AS A WILDCARD
        unmapped_reads = expand( "analysis/STAR/{{sample}}/{{sample}}.Unmapped.out.{mate}", mate=_mates),
        sjtab = "analysis/STAR/{sample}/{sample}.SJ.out.tab",
    params:
        stranded=strand_command,
        gz_support=gz_command,
        prefix=lambda wildcards: "analysis/STAR/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
        keepPairs = _keepPairs,
        star_index = config['star_index'],
    threads: 8
    message: "Running STAR Alignment on {wildcards.sample}"
    priority: 10
    benchmark:
        "benchmarks/{sample}/{sample}.run_STAR.txt"
    shell:
        "STAR --runMode alignReads --runThreadN {threads}"
        " --chimSegmentMin 12"  # added b/c it is specified essential to invoke chimeric read detection & reporting by official STAR-Fusion wiki (https://github.com/STAR-Fusion/STAR-Fusion/wiki)
        " --chimOutJunctionFormat 1"
        " --genomeDir {params.star_index}"
        " --readFilesIn {input} {params.gz_support}" 
        " --outFileNamePrefix {params.prefix}."
        " --outSAMstrandField intronMotif"
        " --outSAMmode Full --outSAMattributes All {params.stranded}"
        " --outSAMattrRGline {params.readgroup}"
        " --outSAMtype BAM SortedByCoordinate"
        " --limitBAMsortRAM 50000000000" #50 Gib
        " --quantMode GeneCounts"
        " --outReadsUnmapped Fastx"
        " --outSAMunmapped Within {params.keepPairs}"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"
        " && mv {params.prefix}.ReadsPerGene.out.tab {output.counts}"

# index sorted bam file (output of STAR) with samtools index
rule index_bam:
    """INDEX the {sample}.sorted.bam file"""
    input:
        "analysis/STAR/{sample}/{sample}.sorted.bam"
    output:
        "analysis/STAR/{sample}/{sample}.sorted.bam.bai"
    message: "Indexing {wildcards.sample}.sorted.bam"
    benchmark:
        "benchmarks/{sample}/{sample}.index_bam.txt"
    shell:
        "samtools index {input}"

rule generate_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
        star_gene_count_files=expand( "analysis/STAR/{sample}/{sample}.counts.tab", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/" + config["token"] + "/STAR/STAR_Align_Report.csv",
        png="analysis/" + config["token"] + "/STAR/STAR_Align_Report.png",
        gene_counts="analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv"
    message: "Generating STAR report"
    #priority: 3
    benchmark:
        "benchmarks/" + config["token"] + "/generate_STAR_report.txt"
    run:
        log_files = " -l ".join( input.star_log_files )
        count_files = " -f ".join( input.star_gene_count_files )
        shell( "perl viper/modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript viper/modules/scripts/map_stats.R {output.csv} {output.png}" )
        shell( "perl viper/modules/scripts/raw_and_fpkm_count_matrix.pl -f {count_files} 1>{output.gene_counts}" )

rule batch_effect_removal_star:
    input:
        starmat = "analysis/" + config["token"] + "/STAR/STAR_Gene_Counts.csv",
        annotFile = config["metasheet"]
    output:
        starcsvoutput="analysis/" + config["token"] + "/STAR/batch_corrected_STAR_Gene_Counts.csv",
        starpdfoutput="analysis/" + config["token"] + "/STAR/star_combat_qc.pdf"
    params:
        batch_column="batch",
        datatype = "star",
        token = config['token'],
    message: "Removing batch effect from STAR Gene Count matrix, if errors, check metasheet for batches, refer to README for specifics"
    #priority: 2
    benchmark:
        "benchmarks/" + config["token"] + "/batch_effect_removal_star.txt"
    shell:
        "Rscript viper/modules/scripts/batch_effect_removal.R {input.starmat} {input.annotFile} "
        "{params.batch_column} {params.datatype} {output.starcsvoutput} {output.starpdfoutput} "
        " && mv {input.starmat} analysis/{params.token}/STAR/without_batch_correction_STAR_Gene_Counts.csv"


rule run_STAR_fusion:
    input:
        bam="analysis/STAR/{sample}/{sample}.sorted.bam", #just to make sure STAR output is available before STAR_Fusion
        left_fastq = getLeftFastq,
        right_fastq = getRightFastq
    params:
        genome_lib_dir = config['genome_lib_dir']
    output:
        #LEN: TODO- CLEAR this out!
        #"analysis/STAR_Fusion/{sample}/star-fusion.fusion_predictions.tsv",
        #"analysis/STAR_Fusion/{sample}/star-fusion.fusion_predictions.abridged.tsv",
        #"analysis/STAR_Fusion/{sample}/FusionInspector-inspect/finspector.FusionInspector.fusions.tsv",
        #"analysis/STAR_Fusion/{sample}/FusionInspector-inspect/finspector.FusionInspector.fusions.abridged.tsv",
        "analysis/STAR_Fusion/{sample}/star-fusion.fusion_predictions.abridged.coding_effect.tsv",
    log:
        "analysis/STAR_Fusion/{sample}/{sample}.star_fusion.log"
    message: "Running STAR fusion on {wildcards.sample}"
    threads: 8
    benchmark:
        "benchmarks/{sample}/{sample}.run_STAR_fusion.txt"
    shell:
        "STAR-Fusion --chimeric_junction analysis/STAR/{wildcards.sample}/{wildcards.sample}.Chimeric.out.junction "
        "--FusionInspector inspect --left_fq {input.left_fastq} --right_fq {input.right_fastq} "
        "--CPU {threads} "
        "--examine_coding_effect "
        "--genome_lib_dir {params.genome_lib_dir} --output_dir analysis/STAR_Fusion/{wildcards.sample} >& {log}"
        #DROPPING THESE- output should be star-fusion.fusion_predictions.tsv
        # and star-fusion.fusion_predictions.abridged.tsv
        
        #THESE names are very mouch tied to the version of STAR-Fusion
        #These are for STAR-Fusion=1.5.0
        #" && cp analysis/STAR_Fusion/{wildcards.sample}/FusionInspector-inspect/finspector.fusion_predictions.final {output[0]}"
        #" && cp analysis/STAR_Fusion/{wildcards.sample}/FusionInspector-inspect/finspector.fusion_predictions.final.abridged {output[1]}"
        #" && touch {output[1]}" # For some sample, final.abridged is created but not .final file; temp hack before further investigate into this


rule run_STAR_fusion_report:
    """NOTE: this rule is broken b/c I'm not sure which output file to use
    from run_STAR_fusion"""
    input:
        #LEN: TODO- CLEAR this out!
        #sf_list = expand("analysis/STAR_Fusion/{sample}/{sample}.fusion_predictions.final.abridged", sample=config["ordered_sample_list"]),
        #sf_list = expand("analysis/STAR_Fusion/{sample}/star-fusion.fusion_predictions.abridged.tsv", sample=config["ordered_sample_list"]),
        #sf_list = expand("analysis/STAR_Fusion/{sample}/FusionInspector-inspect/finspector.FusionInspector.fusions.tsv", sample=config["ordered_sample_list"]),
        sf_list = expand("analysis/STAR_Fusion/{sample}/star-fusion.fusion_predictions.abridged.coding_effect.tsv", sample=config["ordered_sample_list"]),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/" + config["token"] + "/STAR_Fusion/STAR_Fusion_Report.csv",
        png="analysis/" + config["token"] + "/STAR_Fusion/STAR_Fusion_Report.png"
    message: "Generating STAR fusion report"
    benchmark:
        "benchmarks/" + config["token"] + "/run_STAR_fusion_report.txt"
    shell:
        "python viper/modules/scripts/STAR_Fusion_report.py -f {input.sf_list} 1>{output.csv} "
        "&& Rscript viper/modules/scripts/STAR_Fusion_report.R {output.csv} {output.png}"


rule run_rRNA_STAR:
    input:
        getFastq
    output:
        bam=protected("analysis/STAR_rRNA/{sample}/{sample}.sorted.bam"),
        log_file="analysis/STAR_rRNA/{sample}/{sample}.Log.final.out"
    params:
        stranded=rRNA_strand_command,
        gz_support=gz_command,
        prefix=lambda wildcards: "analysis/STAR_rRNA/{sample}/{sample}".format(sample=wildcards.sample),
        readgroup=lambda wildcards: "ID:{sample} PL:illumina LB:{sample} SM:{sample}".format(sample=wildcards.sample),
        star_rRNA_index=config['star_rRNA_index'],
    threads: 8
    message: "Running rRNA STAR for {wildcards.sample}"
    benchmark:
        "benchmarks/{sample}/{sample}.run_rRNA_STAR.txt"
    shell:
        "STAR --runMode alignReads --runThreadN {threads}"
        " --genomeDir {params.star_rRNA_index}"
        " --readFilesIn {input} {params.gz_support}"
        " --outFileNamePrefix {params.prefix}."
        " --outSAMmode Full --outSAMattributes All {params.stranded}"
        " --outSAMattrRGline {params.readgroup}"
        " --outSAMtype BAM SortedByCoordinate"
        " --limitBAMsortRAM 45000000000"
        " && mv {params.prefix}.Aligned.sortedByCoord.out.bam {output.bam}"

rule index_STAR_rRNA_bam:
    """INDEX the STAR_rRNA/{sample}.sorted.bam file"""
    input:
        "analysis/STAR_rRNA/{sample}/{sample}.sorted.bam"
    output:
        "analysis/STAR_rRNA/{sample}/{sample}.sorted.bam.bai"
    message: "Indexing STAR_rRNA {wildcards.sample}.sorted.bam"
    benchmark:
        "benchmarks/{sample}/{sample}.index_STAR_rRNA_bam.txt"
    shell:
        "samtools index {input}"

rule generate_rRNA_STAR_report:
    input:
        star_log_files=expand( "analysis/STAR_rRNA/{sample}/{sample}.Log.final.out", sample=config["ordered_sample_list"] ),
        force_run_upon_meta_change = config['metasheet'],
        force_run_upon_config_change = config['config_file']
    output:
        csv="analysis/" + config["token"] + "/STAR_rRNA/STAR_rRNA_Align_Report.csv",
        png="analysis/" + config["token"] + "/STAR_rRNA/STAR_rRNA_Align_Report.png"
    message: "Generating STAR rRNA report"
    benchmark:
        "benchmarks/" + config["token"] + "/run_rRNA_STAR_report.txt"
    run:
        log_files = " -l ".join( input.star_log_files )
        shell( "perl viper/modules/scripts/STAR_reports.pl -l {log_files} 1>{output.csv}" )
        shell( "Rscript viper/modules/scripts/map_stats_rRNA.R {output.csv} {output.png}" )

rule align_SJtab2JunctionsBed:
    """Convert STAR's SJ.out.tab to (tophat) junctions.bed BED12 format"""
    input:
        "analysis/STAR/{sample}/{sample}.SJ.out.tab"
    output:
        "analysis/STAR/{sample}/{sample}.junctions.bed"
    benchmark:
        "benchmarks/{sample}/{sample}.align_SJtab2JunctionsBed.txt"
    shell:
        "viper/modules/scripts/STAR_SJtab2JunctionsBed.py -f {input} > {output}"
