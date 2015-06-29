#!/usr/bin/env python

import sys
import argparse
import logging
from itertools import islice
from collections import defaultdict

# takes GFF annotation + SAM alignment file and outputs a table of insertions
# (mapped reads) by gene.
#

def parse_args():
    ''' return arguments
        >>> args = parse_args()
    '''

    parser = argparse.ArgumentParser()

    parser.add_argument('--log', default='/dev/stderr',
                        help='log file (default=stderr)')
    parser.add_argument('--gff', help='GFF annotation file')
    parser.add_argument('--sam', help='Aligment in SAM format')
    parser.add_argument('--out', help='output. default = stdout',
                                    default='/dev/stdout')

    return parser.parse_args()


def load_genes(handle):
    ''' parse a GFF file return an iterable of genes and their coordinates:

        [ { 'start': int, 'stop': int, 'desc': 'description' } ]

        GFF uses 1-based indexing

    '''

    for line in handle:
        if line.startswith('#'):
            continue

        line = line.strip().split("\t")

        _type, start, stop, desc = \
                line[2], \
                int(line[3]), \
                int(line[4]), \
                line[-1]

        desc = desc.split(';')[2].lstrip('Name=')

        yield { 'type': _type,
                'start': start,
                'stop': stop,
                'desc': desc }


def load_sam(handle):
    ''' load_sam
    
        SAM uses 1-based indexing 
        Calculates stop position based on sequence length
    '''

    for line in handle:

        if line.startswith('@'):
            continue

        line = line.strip().split("\t")

        start = int(line[3])
        bit = int(line[1])
        strand = { True: '+', False: '-' }[bit == 0x00 and bit != 0x16]
        seq   = line[9]

        if strand == '+':
            stop = start + len(seq)
        elif strand == '-':
            stop = start - len(seq)

        yield { 'start': start,
                'stop': stop,
                'length': len(seq),
                'strand': strand }


def main():
    '''
        >>> main() # stuff happens
    '''

    args = parse_args()
    logging.basicConfig(filename=args.log, level=logging.INFO)

    # load list of genes

    logging.info('loading genes from GFF: %s' % args.gff)

    with open(args.gff) as handle:
        genes = [ g for g in load_genes(handle) if g['type'] == 'gene' ]

    logging.info('loaded %s genes' % len(genes))

    logging.info('loading sam from %s' % args.sam)

    with open(args.sam) as handle:
        inserts = list((i['start'], i['stop']) for i in load_sam(handle))

    logging.info('loaded %s unique inserts' % len(inserts))

    logging.info('this will take a long time... get some coffee')

    # gene -> unique inserts
    gene_to_unique = defaultdict(lambda: defaultdict(int))

    # O(too many)
    for g in genes:
        #logging.info(g['desc'])
        # find start and stop positions of gene, sort in case the gene
        # is on the opposite strand (stop < start)
        a, b = min((g['start'], g['stop'])), max((g['start'], g['stop']))

        # determine if the insert is within the range of the gene. Do not
        # remove inserts after they've been found within a gene as genes can
        # overlap
        for i, j in inserts:
            if (i >= a) and (i <= b) \
                    or \
               (j >= a) and (j <= b):

                if i < j:
                    strand = '+'
                else:
                    strand = '-'
                # only care about start position (different stop positions
                # aren't considered unique inserts because reads are
                # trimmed/different lengths)
                key = '%s,%s' % (i, strand)
                # count number of unique inserts by gene
                # each key should be a unique insert
                gene_to_unique[g['desc']][key] += 1


    # print table of insert location (gene) -> number of reads
    with open(args.out, 'w') as handle:
        print >> handle, 'gene,total_inserts,site,strand,reads'
        for gene, inserts in gene_to_unique.iteritems():
            for insert, reads in inserts.iteritems():
                print '%s,%s,%s,%s' % (gene, len(inserts), insert, reads)

    # need another table of insert location (gene) -> number of unique inserts

if __name__ == '__main__':
    main()
