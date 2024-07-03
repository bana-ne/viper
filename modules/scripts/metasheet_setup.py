#!/usr/bin/env python

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

#-------------------------------------
# @author: Mahesh Vangala
# @email: vangalamaheshh@gmail.com
# @date: July, 1st, 2016
#-------------------------------------
# =============================================================================
# Copyright: Copyright 2024, viper_adaption
# License: GPL
# Version: v0.1
# Editor: Vanessa Mandel
# Email: schmollv@mytum.de
# Created: 29.05.2024 (DD.MM.YYYY)
# Last Modified: 
# =============================================================================
""" Description:
    Update metadata information in Snakemake's config dict
    - remove invalid characters and convert to unix line endings inside the metasheet csv file
    - Extract column names starting with `comp_`
    - add file name(s) for each sample to the config dict
        - edited by Vanessa Mandel: Give option to only specify a path to fastq files (config["fastq_path_prefix"]) in the config.yaml file
            or specify the sample names and samples (see commented option in standard config.yaml file)
        - if only path is specified the SampleName column from the metasheet is used to search for existing files (i.e. config["fastq_path_prefix"]{sample}/{sample}*.fastq(.gz)?)
"""

import pandas as pd
from collections import defaultdict
import warnings
import glob

def updateMeta(config):
    _sanity_checks(config)
    metadata = pd.read_table(config['metasheet'], index_col=0, sep=',', comment='#') # read metadata csv file
    config["comparisons"] = [c[5:] for c in metadata.columns if c.startswith("comp_")] # extract column names starting with `comp_` and remove the five characters `comp_`
    config["comps"] = _get_comp_info(metadata)
    config["metacols"] = [c for c in metadata.columns if c.lower()[:4] != 'comp'] # add names of the metadata columns to config
    # check if fastq_path_prefix is specified in the config file under the samples key
    if "fastq_path_prefix" in config["samples"]:
        _update_sample_file_info(config,metadata)
    else : # assume that sample names and their corresponding files are specified in config file
        config["file_info"] = { sampleName : config["samples"][sampleName] for sampleName in metadata.index } # extract info with file paths specified in the config.yaml file for all samples specified in the metadata csv file
    config["ordered_sample_list"] = metadata.index
    return config


def _sanity_checks(config):
    #metasheet pre-parser: converts dos2unix line endings, catches invalid chars
    _invalid_map = {'\r':'\n', '(':'.', ')':'.', ' ':'_', '/':'.', '$':''}
    _meta_f = open(config['metasheet'])
    _meta = _meta_f.read()
    _meta_f.close()

    _tmp = _meta.replace('\r\n','\n')
    #check other invalids
    for k in _invalid_map.keys():
        if k in _tmp:
            _tmp = _tmp.replace(k, _invalid_map[k])

    #did the contents change?--rewrite the metafile
    if _meta != _tmp:
        #print('converting')
        _meta_f = open(config['metasheet'], 'w')
        _meta_f.write(_tmp)
        _meta_f.close()


def _get_comp_info(meta_info):
    ''' create dictionary with comparison name (respective column name with out `comp_` prefix) and the sample names corresponding to control (1) and treatment (2) groups. 
        comps_info looks something like this: {comp_name:{'control':['sample1','sample3',...], 'treat':['sample2',...]}, comp_name:...}
    '''
    comps_info = defaultdict(dict)
    for comp in meta_info.columns:
        if comp[:5] == 'comp_':
            comps_info[comp[5:]]['control'] = meta_info[meta_info[comp] == 1].index
            comps_info[comp[5:]]['treat'] = meta_info[meta_info[comp] == 2].index
    return comps_info

def _update_sample_file_info(config,metadata):
    """ dynamically locate and use FASTQ files based on the metadata and specified directory prefix and add it to config
    """
    # extract the prefix
    prefix = config["samples"]["fastq_path_prefix"]
    # Initialize samples dictionary if not already present
    if len(config["samples"]) == 1:
        pass
    else:
        warnings.warn("Additional entries in the config file for samples other than `fastq_path_prefix`. Will extract the sample names from the metasheet and search for the respective files specified in the `fastq_path_prefix`.", UserWarning)
    # iterate over samples, extract file names and save them in the config["samples"] dictionary
    config["samples"] = {} # empty all entries in samples
    for sample in metadata.index:
        search_str = "/".join((prefix+sample,f"{sample}*.fastq.*"))
        fastq_files = glob.glob(search_str)
        # check if more than 2 fastq files exist per sample and if so concatenate all R1 and R2 fastq files with comma 
        if len(fastq_files) > 2:
            leftFastq = ",".join([x for x in fastq_files if "R1" in x])  # TODO check if lanes are always ordered correctly??
            rightFastq = ",".join([x for x in fastq_files if "R2" in x])
            fastq_files = [leftFastq,rightFastq]
        config["samples"][sample] = fastq_files
