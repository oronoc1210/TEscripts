import os
from datetime import datetime

outf = open("batch_data.txt", "w+")

for subdir, dirs, files in os.walk("Sbicolor_libraries"):
    for file in files:
        filepath = subdir + os.sep + file
        if filepath.endswith("fastq.gz"):
            timestamp = os.lstat(filepath).st_mtime
            time = datetime.fromtimestamp(timestamp)
            fileparts = file.split(".")
            outbase = fileparts[0]
            outf.write("%s\t%s\t%s\n" % (subdir, file, outbase))
