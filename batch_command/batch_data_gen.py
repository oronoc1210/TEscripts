import os

outf = open("batch_data.txt", "w+")
tepath = "/global/projectb/scratch/cmodonog/TransposableElements"

for dir, subdirs, files in os.walk("%s/Sbicolor_libraries/" % tepath):
    for file in files:
        filepath = dir + os.sep + file
        if "RTx430_Leaf_Irr1_Rep1" in filepath:
            continue
        if filepath.endswith("fastq.gz"):
           fileparts = filepath.split("/")
           filename = fileparts[-1]
           filedir = "/".join(fileparts[:-1])
           libname = fileparts[-2]
           jamo_id = filename.split(".")[0]
           outf.write("{}\t{}\t{}\t{}\n".format(filedir, filename, libname, jamo_id))
