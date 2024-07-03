#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#---------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#---------------------------
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
    Main Snakefile edited from the Viper pipeline for RNA-seq data processing
    Cornwell M, et al. VIPER: Visualization Pipeline for RNA-seq, a Snakemake workflow for efficient and complete RNA-seq analysis. BMC Bioinformatics. 2018 PMID: 29649993.

    This script includes all the modules defined in the Snakefiles inside the rules subdirectory.
    An overview of the available modules can be found in the README.md file.
"""
# =============================================================================

from modules.scripts.config_setup import updateConfig
from modules.scripts.metasheet_setup import updateMeta
from modules.scripts.utils import getTargetInfo

#---------  CONFIG set up  ---------------
configfile: "config.yaml"   # This makes snakemake load up yaml into config 
config = updateConfig(config)
config = updateMeta(config)
#-----------------------------------------

rule target:
    input: getTargetInfo(config), "analysis/" + config["token"] + "/" + config["token"] + ".html"
    message: "Compiling all output"        

include: "./modules/align.snakefile"         # rules specific to STAR and Fusion
include: "./modules/cuff.snakefile"          # cufflinks' and fpkm plot rules
include: "./modules/readQC.snakefile"        # RSeQC and Picard rules 
include: "./modules/file_format.snakefile"   # bam to bigwig
include: "./modules/snp.snakefile"           # variant calling rules
include: "./modules/preprocess.snakefile"    # preprocess rules for cleaning data up 
include: "./modules/cluster.snakefile"       # PCA, Heatmaps (Sample-Sample & Sample-Feature)
include: "./modules/DE.snakefile"            # DESeq2, Limma and volcano plot rules
include: "./modules/pathway.snakefile"       # GO and KEGG rules
include: "./modules/rsem.snakefile"          # RSEM module
include: "./modules/seurat.snakefile"        # seurat - tSNE rules
include: "./modules/virusseq.snakefile"      # Virusseq rules
include: "./modules/immunology.snakefile"    # Immunology module rules
include: "./modules/final_report.snakefile"  # rules for HTML report 
include: "./modules/cdr3.snakefile"          # CDR3 Analysis module
include: "./modules/gsea.snakefile"          # GSEA Analysis module
include: "./modules/optitype.snakefile"      # HLA Analysis module
