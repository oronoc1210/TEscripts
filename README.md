# TEscripts
Scripts used for transposable element expression analysis in Sorghum bicolor.
Work from October 2018.

## Author

Conor O'Donoghue cmodonoghue@lbl.gov

## Requirements

Most scripts require numpy and pandas.

STAR aligner from https://github.com/alexdobin/STAR was used to align sample libraries to the reference genome,
and TEcount from https://github.com/mhammell-laboratory/tetoolkit/ was used to count TE expression.

# script logic

What I had to start with were two excel sheets:
  1. samples_reference.xlsx : contained data about each sample. sample name, tissue type, treatment, date, etc
  2. libraries.txt : A lookup table for what each sample is called within jamo (the JGI's mongodb-based data management system)
  
The Sorghum bicolor reference genome

Two genome annotations:
  1. Sorghum bicolor gene annotation in gff3 format
  2. Sorghum bicolor transposable element annotation in gff3 format
  
and a fasta file of Transposable element sequences, whose headers contained more detailed information
of the TE classifiation than what was present in the gff3 file.

## fetch_all_libraries.py

The JGI Archive and Metadata Organizer (jamo) is the data management system we use to fetch archived data from NERSC tapes.
It requires a library name (or a combination of metadata queries) to find the correct library, which are luckily located 
in libraries.txt

This script used the sample metadata in samples_reference.xlsx and library name from libraries.txt to query jamo using its web API, 
restore the desired fastq, and then create a symbolic link to it in the the correct subdirectory based on sample metadata.

Searching each library in jamo returns a few different files (some processing is done on each library), but this script makes sure to
only fetch and create a link to the raw fastq file.

## gff3_to_gtf.py

tetoolkit has two scripts -- TEtranscripts and TEcount. Both count the expression of transposable elements,
but TEtranscripts takes multiple libraries as arguments and then also performs differential expression analysis
between the libraries given using the count data using DESeq2.

We wanted to do our own analysis, so merely generating the count tables for each library was sufficient.
Hence TEcount is sufficient.

Both scripts, however, require the annotation files -- for both gene and transposable elements -- to be in gtf format, not gff3.

The solution for the gene annotation is easy. The gffread command from the cufflinks suite easily converts files
from gff3 to gtf from the command line.

However, this doesn't work for the transposable element annotation: tetoolkit requires specific headers in the comment field
of the transposable element gtf, which gffread doesn't do. It requires TE gene_id, transcript_id, family_id, and class_id.
Only the name of the transposable element was present in the gff3 file (which is used for gene_id), but the headers in the 
Sorghum bicolor transposable elements fasta file included family and class information. Transcript_id was generated simply
from adding a number into the end of the gene_id, to be able to differentiate between different transcripts of the same TE.

So, gff3_to_gtf.py takes the TE gff3 and fasta files, and writes a new gtf file with the correct headers by using the TE name
in the gff3 file to look up the rest of the data in the fasta file. The result is a gtf file that works with TEcount.

# batch commands

The Berkeley Lab computing cluster uses the SLURM workload manager. For each library I needed to map the fastq to the reference
genome using STAR, and then count transposable element expression using TEcount. This requires ~40G of memory and takes ~1hr each.
So, I needed to set up an array command, which would submit a job for each library.

## batch_data.txt

Each line contains the arguments for one job. The command I submit will iterate over this.

Column1: filepath 
Column2: filename
Column3: output base name

## batch_data_gen.py

A simple command that went into the directory where I kept all of the sample fastq's, 
and walked through all of the subdirectories looking for fastq's. For each file it found,
it would write to batch_data.txt the path to the directory containing the fastq as the filepath, 
the name of the fastq as the filename, and the name of the sample as the output base name.

## batch_wrapper.sh

This is the command that will be done for each job. Every library individaully will be mapped to the reference
using the STAR aligner, and then have transposable element expression counted using TEcount.

## batch_cmd.sl

SLURM file that actually handles submitting the array of jobs. This is the mechanism that iterates over 
batch_data.txt and submits a new run of batch_wrapper.sh for each line.

To kick off the whole process, the batch_cmd.sl file is submitted to SLURM using the sbatch command,
with the specified range of the array used.

For a testrun, I ran sbatch --array=1-8 batch_cmd.sl
with the time for each run set in batch_cmd to 2 hours to start with.
I averaged the time it actually took each job, doubled that time, and then used that as the time each job is given
for the rest of the libraries. I then ran the rest of the libraries with --array=9-809 batch_cmd.sl
