#!/usr/bin/env python
# vim: syntax=python tabstop=4 expandtab

# -*- coding: utf-8 -*-
# =============================================================================
# Author: Vanessa Mandel
# Copyright: Copyright 2024, viper_adaption
# License: GPL
# Version: v0.1
# Maintainer: Vanessa Mandel
# Email: schmollv@mytum.de
# Created: 28.05.2024 (DD.MM.YYYY)
# Last Modified: 
# =============================================================================
""" Description:
    This is a wrapper file that selects the correct rules to run for the alignment, or whether the alignment was already performed 
    with another tool. The purpose is to select the right set of rules by either running the selected alignment tool, or by formatting the files from different tools.
"""
# =============================================================================

#TODO pick the correct alignment snakefile depending on the input files.
# If fastq as input run align_star_fusion.snakefile
# If evadb output run align_evadb.snakefile
include: "./align_star_fusion.snakefile"         # rules specific to STAR and Fusion
#include: "./modules/align_evadb.snakefile"         # rules specific to make the EVAdb output files work with the rest of the pipeline