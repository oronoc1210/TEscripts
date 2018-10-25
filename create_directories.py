import requests
import os
import sys
import logging

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('create_directories.log')
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

logger.info('Starting...')
# Yes, for every single entry. ~1200. Might take a while.
for index, row in df_md.iterrows():
    # Gathering sample data
    sample_name = row['Sample name']
    try:
        sample_lib = df_lib.loc[df_lib["5=sampleName"] == sample_name].iloc[0][0]
    except IndexError:
        logger.error("Sample %s not in LIBRARIES_TABLE.txt, no library name." % sample_name)
        sample_lib = False
    tissue = tissues[row['Tissue']]
    genotype = row['Genotype']
    treatment = row['Treatment'].split()[1]
    replicate = row['Replicate']
    logger.info("Name:%s\tLib:%s\tGenotype:%s\tTissue:%s\tTreatment:Irr%s\tReplicate:%s"
                % (sample_name, sample_lib, genotype, tissue, treatment, replicate))

    # Making sure output path exists and then creating symbolic link to data
    directory = "%s/Sbicolor_libraries/%s_%s_Irr%s_Rep%s" % (te_path, genotype, tissue, treatment, replicate)
    if sample_lib:
        subdir = "%s/%s_%s" % (directory, sample_lib, sample_name)
    else:
        subdir = "%s/%s" % (directory, sample_name)
    if not os.path.exists(directory):
        logger.info("directory %s doesn't exist, creating..." % directory)
        os.mkdir(directory)
    if not os.path.exists(subdir):
        os.mkdir(subdir)
    else:
        logger.info("Directory %s already exists." % subdir)
