import sys
import os
import random
import argparse
import time
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import numpy as np

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
    ids = make_ids(args.i)
    good, bad = cut_runts(ids)
    clusters = cluster(set(good.keys()))
    centers = center(clusters, good)
    if args.graph:
        jackpottogram(clusters, good)
        #bcodePCluster(clusters)
        #minDistHist(centers)
    #print_cids(centers, clusters, good, args.out+".cid")

#reads fasta/q file to get set of barcodes to cluster
#barcodes: {'a', 'b', 'c'}, no repeats
#ids: {'a': ['a1', 'a2'], 'b': ['b1', 'b2', 'b3']}, keys will not be empty
def make_ids(f):
    ids = {}
    i = 0
    start = time.time()
    with open(f, 'r') as f:
        for lab, seq, exp in readfx(f):
            if seq.upper() in ids:
                ids[seq.upper()] += [lab.split()[0]]
                #ids[seq.upper()] += 1
            else:
                ids[seq.upper()] = [lab.split()[0]]
                #ids[seq.upper()] = 1
            i += 1
            if time.time() - start > 1 and args.verbose:
                sys.stdout.write("\rRead %i sequences" %i)
                sys.stdout.flush()
                start = time.time()
        if args.verbose:
            sys.stdout.write("\rRead all %i sequences\tFound %i unique sequences\n" %(i, len(ids)))
    return ids

def cut_runts(ids):
    g = ids.copy()
    b = {}
    k = 0
    c = 0
    start = time.time()
    for barcode in ids:
        if len(ids[barcode]) <= args.cut:
            b[barcode] = g.pop(barcode)
            k -= 1
            c += 1
        k += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rKept: %i\tCut: %i" %(k, c))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rKept: %i\tCut all %i barcodes smaller than %i\n" %(k, c, args.cut))
    return g, b

#takes set of barcodes, outputs array of string array of clusters
#No string array will be empty
#[["a", "b", "c"], ["d"], ["e", f"], ["g", "h", "i", "j"]]
def cluster(barcodes):
    clusters = []
    b = barcodes
    i = 0
    start = time.time()
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
        i += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rGrouped %i clusters" %i)
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rGrouped %i clusters\n" %i)
    return clusters

#takes an array of string array of clusters, and dictionary of unique barcode
#frequencies, returns a string array of "centers" of each cluster
#["c0", "c1", "c2"]
def center(clusters, ids):
    centers = []
    check = set()
    bases = ['A', 'C', 'T', 'G']
    t = 0
    start = time.time()
    for cluster in clusters:
        distr = []
        for barcode in cluster:
            n = len(ids[barcode])
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
        t += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rCentr'd %i clusters" %t)
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rCentered all %i clusters\n" %t)
    return centers

def print_cids(centers, clusters, ids, f):
    with open(f, 'w') as f:
        labs  = []
        t = 1
        start = time.time()
        for i in clusters.pop(0):
            labs += ids[i]
        f.write(">" + centers.pop(0) + "\n" + str(labs))
        while centers and clusters:
            labs  = []
            for i in clusters.pop(0):
                labs += ids[i]
            f.write("\n>" + centers.pop(0) + "\n" + str(labs))
            t += 1
            if time.time() - start > 1 and args.verbose:
                sys.stdout.write("\rWrote %i clusters" %t)
                sys.stdout.flush()
                start = time.time()
        if args.verbose:
            sys.stdout.write("\rWrote all %i clusters\n" %t)
        if centers:
            print "More centers than clusters"
        if clusters:
            print "More clusters than centers"

#takes barcode, outputs set of all 1 off barcodes
#{'a', 'b', 'c'}, no repeats
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

def jackpottogram(clusters, ids):
    n = []
    total = 0
    for cluster in clusters:
        x = 0
        for barcode in cluster:
            x += len(ids[barcode])
        n.append(x)
        total += x
    n = sorted(n, reverse = True)
    even = float(total) / len(clusters)
    plot = plt.figure()

    slices = []
    labels = []
    rest = 0
    other = 0
    o = 0
    for rest in range(len(n)):
        per = 100 * (n[rest] / float(total))
        if per > float(even)/2 or rest < 7:
            slices.append(per)
            labels.append('%.3f%% (%i reads)' %(per, n[rest]))
        elif n[rest] == 1:
            break
        else:
            other += n[rest]
            o += 1
    per = 100 * (other / float(total))
    slices.append(per)
    labels.append('Other %i Clusters:\n%.3f%% (%i reads)' %(o, per, other))
    if rest + 1 != len(n):
        ones = 100 * ((len(n)-rest)/float(total))
        slices.append(ones)
        labels.append("1's: %.3f%%" %ones)

    color = cm.Pastel2(np.linspace(0.,1.,len(slices)))

    handles, text = plt.pie(slices, colors=color, startangle = 90)
    for handle in handles:
        handle.set_edgecolor('white')
        handle.set_linewidth(.05)
    plt.legend(handles, labels, title = "Largest Clusters", loc="upper right", prop={'size':8})
    plt.axis('equal')
    plt.title('Reads per Cluster by Percentage\nCut Barcodes <= %i' %args.cut)
    p = PdfPages(args.out.split('.')[0] + "_jkptogram.pdf")
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
    parser.add_argument("-c", "--cut", help = "cut out barcodes with <= reads than this value",
        type = float, default = 1)
    parser.add_argument("-g", "--graph", help = "creates graphs if specified", 
        action = "store_true", default = False)
    parser.add_argument("-v", "--verbose", help = "output progress information to terminal",
        action = "store_true", default = False)
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()