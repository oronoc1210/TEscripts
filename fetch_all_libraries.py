import requests
import os
import sys
import logging
import re
from datetime import datetime

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

te_path = "/global/projectb/scratch/cmodonog/TransposableElements"
# df = dataframe. md = metadata. lib=libraries.
# The 'samples reference' has all data about the sample except the lib name.
# So also need to open LIBRARIES TABLE to look up lib name using sample name.
df_md = pd.read_excel('%s/EPICON_samples_reference.xlsx' % te_path)
df_lib = pd.read_table("%s/LIBRARIES_TABLE.txt" % te_path)
tissues = {"R": "Root", "L": "Leaf", "P": "P"}

class SkipSample(Exception): pass

def datetime_gen(datestring):
    datearr = re.split("-|:|T", datestring)
    datetup = tuple(map(lambda x: int(x), datearr))
    dt = datetime(*datetup)
    return dt


logger.info("Reading in sample data...")
# Yes, for every single entry. ~1200. Might take a while.
for index, row in df_md.iterrows():
    # Gathering sample data
    sample_name = row['Sample name']
    try:
        sample_lib = df_lib.loc[df_lib["5=sampleName"] == sample_name].iloc[0][0]
    except IndexError:
        logger.info("Sample %s not in LIBRARIES_TABLE.txt, no library name." % sample_name)
        sample_lib = False
    tissue = tissues[row['Tissue']]
    genotype = row['Genotype']
    treatment = row['Treatment'].split()[1]
    replicate = row['Replicate']
    logger.info("Name:%s\tLib:%s\tGenotype:%s\tTissue:%s\tTreatment:Irr%s\tReplicate:%s"
                % (sample_name, sample_lib, genotype, tissue, treatment, replicate))
    directory = "%s/Sbicolor_libraries/%s_%s_Irr%s_Rep%s" % (te_path, genotype, tissue, treatment, replicate)
    if sample_lib:
        subdir = "%s/%s_%s" % (directory, sample_lib, sample_name)
    else:
        subdir = "%s/%s" % (directory, sample_name)
    try:
        for dir, sub_dirs, files in os.walk(subdir):
            for file in files:
                if file.endswith("fastq.gz"):
                    filepath = dir + os.sep + file
                    logger.info("Link already present at %s. Skipping..." % (filepath))
                    raise SkipSample()
    except SkipSample:
        continue


    # Finding jamo id
    headers = {"Content-Type": "application/json"}
    jamo_url = "https://sdm2.jgi-psf.org/api/metadata/pagequery"
    query = ("metadata.sow_segment.sample_name=%s" % sample_name)
    data = {"query": query,
            "fields": ["file_path", "file_name", "file_date"]
           }
    r = requests.post(jamo_url, data=data)

    # Grabbing relevant jamo metadata and fetching fastq
    # Goal: getting the MOST RECENT, FILTERED fastq.
    # So, if multiple filtered, will take most recent.
    # If none filtered, will take most recent .fastq file
    json = r.json()
    _id = False
    _time = False
    _name = False
    for item in json['records']:
        # Looking for fastq's. Not chaff.tar's etc
        if not item['file_name'].endswith('fastq.gz'):
            continue
        # If a filter-RNA has already been found, only filter-RNAs can replace it.
        if _name and 'filter-RNA' in _name:
            if 'filter-RNA' not in item['file_name']:
                continue
        # If filter-RNA hasn't been found yet, don't want FAIL or BISULFITE fastqs.
        if 'FAIL' in item['file_name'] or 'BISULFITE' in item['file_name']:
            continue
        # If we already have a record, will replace with this record if it's newer.
        if _time:
            newtime = datetime_gen(item['file_date'])
            oldtime = datetime_gen(_time)
            youngest = max(newtime, oldtime)
            if youngest == oldtime:
                continue
        _id = item['_id']
        _time = item['file_date']
        _path = item['file_path']
        _name = item['file_name']
        fetch_arg = "{files: [%s]}" % _id
    if not _id:     
        logger.error("No file found for this sample! Records as follows: %s"
                     % (str(json['records'])))
        continue
    r2 = requests.post(jamo_url, fetch_arg)

    # Making sure output path exists and then creating symbolic link to data
    fastq_path = "%s/%s" % (_path, _name)
    if sample_lib:
        fastq_link = "%s/%s.filter-RNA.fastq.gz" % (subdir, sample_lib)
    else:
        fastq_link = "%s/%s.filter-RNA.fastq.gz" % (subdir, sample_name)
    if not os.path.islink(fastq_link):
        logger.info("Creating link of %s at %s" % (fastq_path, fastq_link))
        os.symlink(fastq_path, fastq_link)
    else:
        logger.info("Link %s already exists" % (fastq_link))
