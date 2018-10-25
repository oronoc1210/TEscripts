# TEscripts
Scripts used for transposable element analysis.

## author

Conor O'Donoghue cmodonoghue@lbl.gov

## requirements

Most scripts require numpy and pandas.

STAR aligner from https://github.com/alexdobin/STAR was used to align sample libraries to the reference genome,
and TEcount from https://github.com/mhammell-laboratory/tetoolkit/ was used to count TE expression.

# script logic

These scripts were used to perform analysis on transposable element expression in Sorghum Bicolor.
What I had to start with were two excel sheets:
  1. samples_reference.xlsx : contained data about each sample. sample name, tissue type, treatment, date, etc
  2. libraries.xlsx : A lookup table for what each sample is called within jamo (the JGI's mongodb-based data management system)

## create_directories.py

There were about 1200 samples. While not necessary, I wanted to organize them by type.
I created a directory for each combination of genotype, tissue type, treatment, and replicate.
Within each of these directories were subdirectories for each sample.
I organized things this way as one of the main ways I wanted to analyze the results was TE expression over time.
Thus within each type-combination directory, 
there should only be libraries of the same genotype, tissue type, treatment, and replicate but with different sample dates.
