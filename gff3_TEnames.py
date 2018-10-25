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
    dscr="""Change a Transposable Element .gff3 file so that the name
         column is replaced with the 'name:' field in comment column,
         for more informative fasta headers when using bedtools getfasta.
         if no outbase provided, will use same base name as gff3 file.
         """
    ex="""Example command: gff3_to_fa.py --gff example.repeatmasker.gff3
       --filter --gzipout --outbase somethingelse.repeatmasker
       """
    parser = argparse.ArgumentParser(description=dscr, epilog=ex)
    parser.add_argument('--gff', metavar='[TE gff3 path]', dest='gffpath',
                        required=True, help='path to TE gff3 file')
    parser.add_argument('--filter', action='store_true', default=False,
                        dest='filter', help='Filter simple repeats')
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
    if not os.path.isfile(args.gffpath):
        logger.error('%s file does not exist!\n' % (args.gffpath))
        sys.exit(1)
    gff_exts = ('.gff3', '.gff3.gz')
    if not args.gffpath.endswith(gff_exts):
        logger.error('file %s must be in .gff3 format!\n' % (args.gffpath))
        parser.print_help()
        sys.exit(1)
    if not args.outbase:
        for ext in gff_exts:
            if args.gffpath.endswith(ext):
                args.outbase = args.gffpath[:-len(ext)]
    logger.info("gff3 input: %s" % (args.gffpath))
    logger.info("outbase: %s" % (args.outbase))
    if args.gzipout:
        logger.info('Output will be in .gff3.gz format')
    return args

def write_line(line, gzip, outf):
    if gzip:
        outf.write(line.encode())
    else:
        outf.write(line)

def main():
    args = read_args(parser_gen())
    # open gff3 file
    if args.gffpath.endswith('.gz'):
        gff = gzip.open(args.gffpath, mode='rt')
    else:
        gff = open(args.gffpath, mode='r')
    outbase = '_'.join((args.outbase, 'TEnames'))
    if args.gzipout:
        new_gff = gzip.open(''.join((outbase, 'gff3', 'gz')), 'wb')
    else:
        new_gff = open('.'.join((outbase, 'gff3')), 'w')
    gff_lines=1
    logger.info('Changing name column to "name" header...')
    for line in gff:
        # ignore gff3 comments and simple repeats
        if line.startswith('##'):
            write_line(line, args.gzipout, new_gff)
            continue
        if 'Motif:' in line:
            continue
        if args.filter == True:
            no_good = [')n', '-rich']
            if any(x in line for x in no_good):
                continue
        gff_parts = line.split("\t")
        try:
            te_name = gff_parts[8].split(';')[1][5:]
        except:
            logger.error(line)
            sys.exit(1)
        new_line = '\t'.join([
            gff_parts[0],
            gff_parts[1],
            te_name,
            gff_parts[3],
            gff_parts[4],
            gff_parts[5],
            gff_parts[6],
            gff_parts[7],
            gff_parts[8]])
        write_line(new_line, args.gzipout, new_gff)
        if gff_lines % 100000 == 0:
            logger.info('%d features processed' % (gff_lines))
        gff_lines+=1


if __name__ == '__main__':
    main()

