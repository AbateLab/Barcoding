import sys
import os
import random
import argparse
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
    barcodes = None
    counts = None
    start = time.time()
    barcodes, counts = make_set(args.i)
    print "read to set:\t" + str(time.time() - start)
    start = time.time()
    clusters = cluster(barcodes)
    print "cluster:\t" + str(time.time() - start)
    start = time.time()
    centers = center(clusters, counts)
    print "center: \t" + str(time.time() - start)
    start = time.time()
    print_clusters(centers, clusters, args.out+".ctr")
    print "print cl's':\t"+ str(time.time() - start)
    start = time.time()
    if args.rchist:
        readsPCluster(clusters, counts)
        print "readsper:\t" + str(time.time() - start)
        start = time.time()
    if args.bchist:
        bcodePCluster(clusters)
        print "bcodesper:\t" + str(time.time() - start)
        start = time.time()
    if args.mindist:
        minDistHist(centers)
        print "mindist:\t" + str(time.time() - start)
        start = time.time()

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

#takes an array of string array of clusters, and dictionary of unique barcode
#frequencies, returns a string array of "centers" of each cluster
def center(clusters, counts):
    centers = []
    check = set()
    bases = ['A', 'C', 'T', 'G']
    for cluster in clusters:
        distr = []
        for barcode in cluster:
            n = counts[barcode]
            for i in range(len(barcode)):
                while len(distr) <= i:
                    distr.append([0, 0, 0, 0])
                for j in range(len(bases)):
                    if barcode[i] == bases[j]:
                        distr[i][j] += n
        center = ""
        for elem in distr:
            m = max(elem)
            indeces = [i for i, j in enumerate(elem) if j == m]
            center += bases[random.choice(indeces)]
        if center in check:
            print "Multiple Clusters centered at " + center
        check |= {center}
        centers.append(center)
    return centers

def print_cseq(centers, clusters, f):
    with open(f, 'w') as f:
        for center, cluster in zip(centers, clusters):
            f.write(">" + center + "\n" + str(cluster) + "\n")

"""#takes clusters, and centers, and appends center at the end of all labels of
#barcodes in same cluster
def label(clusters, centers):
    #turn cluster's arrays into sets for faster membership checking
    clusters = [set(i) for i in clusters]
    with open(args.i, 'r') as f, open(args.out, 'w') as out:
        for line in f:
            l = line
            seq = f.next()
            if args.fastq:
                f.next()
                f.next()
            if l:
                #assign center to seq
                tag = None
                for i in range(len(clusters)):
                    if seq.rstrip("\n") in clusters[i]:
                        tag = centers[i]
                        break
                out.write(l[:-1] + " " + tag + "\n")
                out.write(seq)"""

#reads fasta/q file to get set of barcodes to cluster
def make_set(f):
    barcodes = set()
    counts = {}
    with open(f, 'r') as f:
        for lab, seq, exp in readfx(f):
            if seq in counts:
                counts[seq] += 1
            else:
                counts[seq] = 1
            barcodes |= {seq.upper()}
    return barcodes, counts

#takes barcode, outputs set of all 1 off barcodes
def neighbors(barcode):
    barcodes = []
    for i in range(len(barcode)):
        for char in ['A', 'C', 'T', 'G']:
            if char != barcode[i]:
                tmp = list(barcode.upper())
                tmp[i] = char
                barcodes.append("".join(tmp))
    return barcodes

#graphs min distances to another "center" for all centers
def minDistHist(centers):
    d = []
    m = []
    for i in range(len(centers)):
        x = []
        tmp = sys.maxint
        for j in range(len(centers)):
            if i == j:
                x.append(sys.maxint)
            elif i < j:
                x.append(sum(c1 != c2 for c1, c2 in zip(centers[i], centers[j])))
            else:
                x.append(d[j][i])
            if x[j] < tmp:
                tmp = x[j]
        d.append(x)
        m.append(tmp)
    p = PdfPages(args.out.split('.')[0] + "_minDist.pdf")
    plot = plt.figure()
    plt.title("Minimum Distances Between Centers")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.xlabel("Minimum Hamming Distance to Another Center")
    plt.ylabel("Number of Centers")
    plt.hist(m, len(centers[0])/3, color = 'blue', alpha = .6)
    plt.xlim(0, len(centers[0]))
    p.savefig(plot)
    p.close()

def readsPCluster(clusters, counts):
    n = []
    for cluster in clusters:
        x = 0
        for barcode in cluster:
            x += counts[barcode]
        n.append(x)
    p = PdfPages(args.out.split('.')[0] + "_RCHist.pdf")
    plot = plt.figure()
    plt.title("Reads per Cluster")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.xlabel("Reads")
    plt.ylabel("Clusters")
    plt.hist(n, 12543, color = 'blue', alpha = .6, cumulative = True, histtype = "step")
    plt.xlim([0, 100000])
    p.savefig(plot)
    p.close()    

def bcodePCluster(clusters):
    p = PdfPages(args.out.split('.')[0] + "_BCHist.pdf")
    plot = plt.figure()
    plt.title("Barcodes per Cluster")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.xlabel("Barcodes")
    plt.ylabel("Clusters")
    plt.xlim([0, 100000])
    plt.hist([len(cluster) for cluster in clusters], 257, color = 'blue', alpha = .6)
    p.savefig(plot)
    p.close()

#ripped from https://github.com/lh3/readfq
def readfx(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="outputs a fasta file of input with group identity\
        appended to labels, reads in fasta file must be same length, and only reads composed solely\
        of A, C, T, G will be considered.")
    parser.add_argument("i", help = "fasta/q file containing sequences to be clustered")
    parser.add_argument("out", help = "output tag")
    parser.add_argument("-md", "--mindist", help = "create hist of minimum hamming distances",
        action = "store_true")
    parser.add_argument("-rc", "--rchist", help = "create hist of reads per cluster",
        action = "store_true")
    parser.add_argument("-bc", "--bchist", help = "create hist of barcodes per cluster",
        action = "store_true")
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()