import sys
import os
import math
import argparse

"""dfsCluster is a barcode clustering program designed to counteract noise
from PCR errors which diffuses each original DNA tag into a cluster of
similar sequences. dfsCluster treats each cluster as a connected component
on a graph of Hamming Space, and connected components can easily be identified
using a depth first search. In treating each cluster as a connected component,
we assume:
1. PCR errors are single base errors
2. Clusters don't diffuse into another cluster (suspected true, not confirmed)
"""

def main():
    if args.fasta:
        print cluster(fa_set(args.input))
    elif args.fastq:
        print cluster(fq_set(args.input))

#takes set of barcodes, outputs array of string array of clusters
def cluster(barcodes):
    clusters = []
    b = barcodes
    while b:
        tmp = b.pop()
        s = [tmp]
        c = [tmp]
        while s:
            for x in neighbors(s.pop()):
                if x in b:
                    s.append(x)
                    c.append(x)
                    b.remove(x)
        clusters.append(c)
    return clusters

#reads fasta file to get set of barcodes to cluster
def fa_set(f):
    barcodes = set()
    with open(f, 'r') as f:
        barcode = ""
        for line in f:
            if line[0] == '>':
                if barcode:
                    barcodes |= {barcode}
                barcode = ""
            else:
                barcode += line.rstrip('\n')
    return barcodes | {barcode}

#reads fastq file to get set of barcodes to cluster
def fq_set(f):
    barcodes = set()
    with open(f, 'r') as f:
        barcode = ""
        go = False
        for line in f:
            if go:
                barcode += line.rstrip('\n')
            if line[0] == '@':
                go = True
            if line[0] == '+' and barcode:
                barcodes |= barcode
                go = False
                barcode = ""
    return barcodes

#takes barcode, outputs set of all 1 off barcodes
def neighbors(barcode):
    barcodes = []
    for i in range(len(barcode)):
        for char in ['a', 'c', 't', 'g']:
            if char != barcode[i]:
                tmp = list(barcode.lower())
                tmp[i] = char
                barcodes.append("".join(tmp))
    return barcodes


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-fa", "--fasta", help = "specify input file as fasta formatted", action="store_true")
    group.add_argument("-fq", "--fastq", help = "specify input file as fastq formatted", action="store_true")
    parser.add_argument("input", help = "input file of barcodes")
    args = parser.parse_args()
    main()