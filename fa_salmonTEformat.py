import sys
import os
import argparse
import logging
import gzip
import re

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

def parser_gen():
    dscr="""Turn a Transposable Element .gff3 file to .gtf format.
         If repeat lib fasta provided, more information can be included in
         gtf header. Otherwise, 'target' header in gff3 file gets
         entered into all gtf headers. Supports I/O in gzip format."""
    ex="""Example command: gff3_to_gtf.py --gff example.repeatmasker.gff3
       --repeatlib example.MIPS_n_Repbase.fa --removesimple"""
    parser = argparse.ArgumentParser(description=dscr, epilog=ex)
    parser.add_argument('--fa', metavar='[TE fa path]', dest='fapath',
                        required=True, help='path to TE fa file')
    parser.add_argument('--repeatlib', metavar='[repeat lib path]',
                        dest='rlibpath', help='path to repeat lib fasta file')
    parser.add_argument('--outbase', metavar='[output base name]',
                        dest='outbase', help="""base name for output file.
                        If no argument given, uses the same base name 
                        as gff3 input.""")
    parser.add_argument('--gzipout', dest='gzipout',
                        action='store_true', default=False,
                        help='if set, output will be in gzip format')
    return parser


def read_args(parser):
    args=parser.parse_args()
    if not os.path.isfile(args.fapath):
        logger.error('%s file does not exist!\n' % (args.gffpath))
        sys.exit(1)
    fa_exts = ('.fa', '.fa.gz')
    if not args.fapath.endswith(fa_exts):
        logger.error('file %s must be in .fa format!\n' % (args.gffpath))
        parser.print_help()
        sys.exit(1)
    if not os.path.isfile(args.rlibpath):
        logger.error('%s file does not exist!\n' % (args.rlibpath))
        sys.exit(1)
    if not args.rlibpath.endswith(('.fa', '.fa.gz')):
        logger.error('file %s must be in .fa format!\n' % (args.rlibpath))
        parser.print_help()
        sys.exit(1)
    if not args.outbase:
        for ext in fa_exts:
            if args.fapath.endswith(ext):
                args.outbase = args.fapath[:-len(ext)]
    logger.info("fa input: %s" % (args.fapath))
    logger.info("repeat lib input: %s" % (args.rlibpath))
    logger.info("outbase: %s" % (args.outbase))
    if args.gzipout:
        logger.info('Output will be in .gtf.gz format')
    return args

def create_lookup(rlibpath, outbase):
    """ If repeat library fasta provided, this parses that data
    and creates a lookup table. """
    logger.info('Generating lookup table from repeat library...')
    if rlibpath.endswith('.gz'):
        rlib = gzip.open(rlibpath, mode='rt')
    else:
        rlib = open(rlibpath, mode='r')
    lookup = dict()
    linecount = 0
    for line in rlib:
        if line.startswith('>'):
            """ headers in this file follow this format:
            >[gene_id]#[family_id]/[class_id] RepbaseID: [gene_id]
            gene_id and family_id are always present, but class_id and
            RepbaseID are not always present.

            example lines in Sorghum:
            >ALFARE1_I#LTR
            >ABR1#DNA\hAT
            >DNA-9-2_Sbi#DNA RepbaseID: DNA-9-2_SBi
            >LINE1-9_SBi#LINE/L1 RepbaseID: LINE1-9_SBi

            Using the gene_id (the only info present in the gff3 file)
            I want to obtain the family_id and class_id
            I split by multiple delimiters using re.split
            and filter out the empty strings from the resulting array.
            filter() returns an object, which I need to turn into a list.

            This gives the following:
            [gene_id, family_id, class_id, 'RepbaseID:', gene_id]
            or a shorter array, if either class_id or RepbaseID are missing.
            (e.g. [gene_id, family_id])
            I don't need the 'RepbaseID:' or the duplicate of gene_id,
            So I remove those.
            They're always the last two elements of the array,
            making removal easy.
            In the case of class_id missing, I put gene_id in its spot
            resulting in [gene_id, family_id, gene_id].
            """
            linecount+=1
            lineparts = list(filter(None, re.split("[/ #'\n']+", line[1:])))
            if 'RepbaseID:' in lineparts:
                lineparts = lineparts[0:-2]
            if '/' not in line:
                lineparts.append(lineparts[0])
            lookup[lineparts[0]] = lineparts[1:]
    logger.info('Lookup table complete.')
    return lookup

def lookup_header(lookup, fa_header, linecount):
    """ If lookup table exists, use this to generate headers
    for the last column in the gtf file."""
    fa_header = fa_header.split()[0]
    te_gene = fa_header[1:]
    if te_gene.endswith(('_I-int', 'int-int')):
        te_gene = te_gene[:-4]
    try:
        te_fam = lookup[te_gene][0]
    except KeyError:
        logger.error('KeyValue error! %s' % (te_gene))
        return 'KeyError'
    te_class = lookup[te_gene][1]
    fasta_header = ">%s_dup%d\t%s\n" % (te_gene, linecount, te_gene)
    return fasta_header

def write_file(line, gzip_out, outf):
    if gzip_out:
        outf.write(line.encode())
    else:
        outf.write(line)


def main():
    args = read_args(parser_gen())
    # open gff3 file
    if args.fapath.endswith('.gz'):
        fa = gzip.open(args.fapath, mode='rt')
    else:
        fa = open(args.fapath, mode='r')
    # open repeat library fasta file, then create lookup dict
    lookup = create_lookup(args.rlibpath, args.outbase)
    # create output file
    outbase = "_".join((args.outbase, "SalmonTEformat"))
    if args.gzipout:
        newfa = gzip.open('.'.join((outbase, 'fa', 'gz')), 'wb')
    else:
        newfa = open('.'.join((outbase, 'fa')), 'w')
    linecount=0
    logger.info('Beginning conversion of fasta headers...')
    for line in fa:
        if not line.startswith(">"):
            write_file(line, args.gzipout, newfa)
            continue
        linecount+=1
        new_header = lookup_header(lookup, line, linecount)
        if new_header == 'KeyError':
            continue
        write_file(new_header, args.gzipout, newfa)
        if linecount % 100000 == 0:
            logger.info('%d features processed' % (linecount))


if __name__ == '__main__':
    main()

