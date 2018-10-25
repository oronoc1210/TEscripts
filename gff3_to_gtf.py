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
    parser.add_argument('--gff', metavar='[TE gff3 path]', dest='gffpath',
                        required=True, help='path to TE gff3 file')
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
    if not os.path.isfile(args.gffpath):
        logger.error('%s file does not exist!\n' % (args.gffpath))
        sys.exit(1)
    gff_exts = ('.gff3', '.gff3.gz')
    if not args.gffpath.endswith(gff_exts):
        logger.error('file %s must be in .gff3 format!\n' % (args.gffpath))
        parser.print_help()
        sys.exit(1)
    if args.rlibpath:
        if not os.path.isfile(args.rlibpath):
            logger.error('%s file does not exist!\n' % (args.rlibpath))
            sys.exit(1)
        if not args.rlibpath.endswith(('.fa', '.fa.gz')):
            logger.error('file %s must be in .fa format!\n' % (args.rlibpath))
            parser.print_help()
            sys.exit(1)
    if not args.outbase:
        for ext in gff_exts:
            if args.gffpath.endswith(ext):
                args.outbase = args.gffpath[:-len(ext)]
    logger.info("gff3 input: %s" % (args.gffpath))
    if args.rlibpath:
        logger.info("repeat lib input: %s" % (args.rlibpath))
    else:
        logger.warn("no repeat lib provided - all headers will be the same!")
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
            lineparts = list(filter(None, re.split("[/ #'\n']+", line[1:])))
            if 'RepbaseID:' in lineparts:
                lineparts = lineparts[0:-2]
            if '/' not in line:
                lineparts.append(lineparts[0])
            lookup[lineparts[0]] = lineparts[1:]
    logger.info('Lookup table complete.')
    return lookup

def lookup_headers(lookup, gff_headers, gff_num):
    """ If lookup table exists, use this to generate headers
    for the last column in the gtf file."""
    te_gene = gff_headers.split(';')[1][5:]
    if te_gene.endswith(('_I-int', 'int-int')):
        te_gene = te_gene[:-4]
    te_tran = '_'.join((te_gene, 'dup%s' % (str(gff_num))))
    try:
        te_fam = lookup[te_gene][0]
    except KeyError:
        logger.error('KeyValue error! %s' % (te_gene))
        return 'KeyError'
    te_class = lookup[te_gene][1]
    gtf_headers = ' '.join([
                  'gene_id "%s";' % (te_gene),
                  'transcript_id "%s";' % (te_tran),
                  'family_id "%s";' % (te_fam),
                  'class_id "%s";\n' % (te_class)])
    return gtf_headers

def main():
    args = read_args(parser_gen())
    # open gff3 file
    if args.gffpath.endswith('.gz'):
        gff = gzip.open(args.gffpath, mode='rt')
    else:
        gff = open(args.gffpath, mode='r')
    # open repeat library fasta file, then create lookup dict
    if args.rlibpath:
        lookup = create_lookup(args.rlibpath, args.outbase)
    # create output file
    if args.gzipout:
        gtf = gzip.open('.'.join((args.outbase, 'gtf', 'gz')), 'wb')
    else:
        gtf = open('.'.join((args.outbase, 'gtf')), 'w')
    gff_lines=1
    logger.info('Beginning conversion from gff3 to gtf...')
    for line in gff:
        # ignore gff3 comments and simple repeats
        no_good = [')n', '-rich', 'Motif:']
        if line.startswith('##') or any(x in line for x in no_good):
            continue
        gff_parts = line.split()
        if args.rlibpath:
            gtf_headers = lookup_headers(lookup, gff_parts[8], gff_lines)
        else:
            te_tran = '_'.join((te_gene, 'dup%s' % (str(gff_lines))))
            gtf_headers = ' '.join([
                          'gene_id "%s";' % (te_gene),
                          'transcript_id "%s";' % (te_tran),
                          'family_id "%s";' % (te_gene),
                          'class_id "%s";\n' % (te_gene)])
        if gtf_headers == 'KeyError':
            continue
        gtf_line = '\t'.join([
                   gff_parts[0],
                   gff_parts[1],
                   "exon",
                   gff_parts[3],
                   gff_parts[4],
                   gff_parts[5],
                   gff_parts[6],
                   gff_parts[7],
                   gtf_headers])
        if args.gzipout:
            gtf.write(gtf_line.encode())
        else:
            gtf.write(gtf_line)
        if gff_lines % 100000 == 0:
            logger.info('%d features processed' % (gff_lines))
        gff_lines+=1


if __name__ == '__main__':
    main()

