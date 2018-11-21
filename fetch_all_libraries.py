import os
import sys
import logging
import re
import subprocess

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

class SkipSample(Exception): 
    pass

class cd:
    """Context manager to temporarily change current working directory"""
    def __init__(self, new_path):
        self.new_path = new_path

    def __enter__(self):
        self.old_path = os.getcwd()
        os.chdir(self.new_path)

    def __exit__(self, etype, evalue, traceback):
        os.chdir(self.old_path)

def datetime_convert(datestring):
    datearr = re.split("-|:|T", datestring)
    datetup = tuple(map(lambda x: int(x), datearr))
    dt = datetime(*datetup)
    return dt

te_path = "/global/projectb/scratch/cmodonog/TransposableElements"
# df = dataframe. md = metadata. lib=libraries.
# The 'samples reference' has all data about the sample except the lib name.
# So also need to open LIBRARIES TABLE to look up lib name using sample name.
with cd(te_path):
    df_md = pd.read_excel('EPICON_samples_reference.xlsx')
    df_lib = pd.read_table('LIBRARIES_TABLE.txt')
tissues = {"R": "Root", "L": "Leaf"}

logger.info("Reading in sample data...")
# ~400 samples. Might take a while.
for index, row in df_lib.iterrows():
    # Gathering sample data
    library_name = row['1=libraryName']
    sample_name = row['5=sampleName']
    """
    tissue = tissues[row['Tissue']]
    genotype = row['Genotype']
    treatment = row['Treatment'].split()[1]
    replicate = row['Replicate']
    """
    epicon_entry = df_md.loc[df_md["Sample name"] == sample_name].iloc[0]
    tissue = tissues[epicon_entry["Tissue"]]
    genotype = epicon_entry["Genotype"]
    treatment = epicon_entry["Treatment"].split()[1]
    replicate = epicon_entry["Replicate"]

    logger.info("Name:%s\tLib:%s\tGenotype:%s\tTissue:%s\tTreatment:Irr%s\tReplicate:%s"
                % (sample_name, library_name, genotype, tissue, treatment, replicate))
    treatment_subdir = "%s/Sbicolor_libraries/%s_%s_Irr%s_Rep%s" % (te_path, genotype, tissue, treatment, replicate)
    library_subdir = "%s/%s" % (treatment_subdir, library_name)

    if not os.path.isdir(treatment_subdir):
        os.mkdir(treatment_subdir)
    if not os.path.isdir(library_subdir):
        os.mkdir(library_subdir)

    # if a fastq already exists in this file location, can move on.
    try:
        for dir, subdirs, files in os.walk(library_subdir):
            for file in files:
                if file.endswith("fastq.gz"):
                    raise SkipSample()
    except SkipSample:
        continue

    jamo_info_cmd = ["jamo", "info", "library", library_name]
    try:
        info_output = subprocess.check_output(jamo_info_cmd)
        logger.info(info_output)
        if "No matching records" in info_output:
            logger.info("Sample has no records in jamo! Skipping...")
            raise SkipSample()
    except subprocess.CalledProcessError as info_err:
        logger.error("JAMO ERROR %s %s" % (info_err.returncode, info_err.output))
        continue
    except SkipSample:
        continue

    # want to create a list of the entries (first [:-1])
    # and then from the last entry (second [-1])
    # take the last piece of data, which is the jamo id (third [-1])
    jamo_entry = info_output.split('\n')[:-1][-1].split()
    jamo_id = jamo_entry[-1]
    jamo_name = jamo_entry[1]
    logger.info("file chosen: jamo name:%s\tjamo id:%s" % (jamo_name, jamo_id))
    # write the jamo metadata in the same directory the link is going.
    # Might need it later, better to grab it now.
    with cd(library_subdir):
        jamo_show_cmd = "jamo show id %s > %s.json" % (jamo_id, library_name)
        jamo_show_cmd = jamo_show_cmd.split()
        try:
            with open("%s.json" % library_name, "w+") as json_file:
                show_output = subprocess.call(jamo_show_cmd, stdout=json_file)
        except subprocess.CalledProcessError as show_err:
            logger.error("JAMO ERROR %s %s" (fetch_err.returncode, fetch_err.output))

    jamo_fetch_cmd = ["jamo", "fetch", "id", jamo_id]
    try:
        fetch_output = subprocess.check_output(jamo_fetch_cmd)
        logger.info(fetch_output)
    except subprocess.CalledProcessError as fetch_err:
        logger.error("JAMO ERROR %s %s" % (fetch_err.returncode, fetch_err.output))
        continue

    # Creating symbolic link to data in correct subdirectory
    with cd(library_subdir):
        logger.info("Creating link of %s at %s" % (jamo_name, library_subdir))
        jamo_link_cmd = ["jamo", "link", "id", jamo_id]
        try:
            link_output = subprocess.check_output(jamo_link_cmd)
            logger.info(link_output)
        except subprocess.CalledProcessError as link_err:
            logger.error("JAMO ERROR %s %s" % (link_err.returncode, link_err.output))
            continue
logger.info("All libraries completed.")
