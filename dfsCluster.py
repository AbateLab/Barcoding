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
    if args.cutcollisions:
        clusters = cutcol(clusters, good)
    centers = center(clusters, good)
    cids = make_cids(clusters, centers, good)
    if args.graph:
        jackpottogram(cids)
        rpc(cids)
        minDistHist(centers)
    print_cids(cids, args.out+".cid")

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
                sys.stdout.write("\rRead {:,} sequences".format(i))
                sys.stdout.flush()
                start = time.time()
        if args.verbose:
            sys.stdout.write("\rRead all {:,} sequences\tFound {:,} unique sequences\n"
                .format(i, len(ids)))
    return ids

def cut_runts(ids):
    g = ids.copy()
    b = {}
    k = 0
    c = 0
    cr = 0
    start = time.time()
    for barcode in ids:
        if len(ids[barcode]) <= args.cutreads:
            b[barcode] = g.pop(barcode)
            k -= 1
            c += 1
            cr += len(ids[barcode])
        k += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rKept: {:,}\tCut: {:,}".format(k, c))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rKept {:,} barcodes, Cut all {:,} barcodes smaller than {:,} ({:,} sequences)\n"
            .format(k, c, args.cutreads, cr))
    return g, b

#takes set of barcodes, outputs array of string array of clusters
#No string array will be empty
#[["a", "b", "c"], ["d"], ["e", f"], ["g", "h", "i", "j"]]
def cluster(barcodes):
    clusters = []
    b = barcodes.copy()
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
            sys.stdout.write("\rGrouped {:,} clusters".format(i))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rGrouped {:,} clusters\n".format(i))
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
            sys.stdout.write("\rCentr'd {:,} clusters".format(t))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rCentered all {:,} clusters\n".format(t))
    return centers

def cutcol(clusters, ids):
    singles = []
    for cluster in clusters:
        top, sec, tot = 0, 0, 0
        for barcode in cluster:
            idc = len(ids[barcode])
            tot += idc
            if idc > top:
                sec, top = top, idc
            elif idc > sec:
                sec = idc
        if float(top - sec) / tot < .72:
            singles.append(cluster)
    return singles

def make_cids(clusters, centers, ids):
    cids = {}
    i = 0
    c = 0
    start = time.time()
    for cn, cl in zip(centers, clusters):
        cids[cn] = []
        for barcode in cl:
            cids[cn] += ids[barcode]
            i += len(ids[barcode])
        c += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rGrouped {:,} ids into {:,} centers".format(i, c))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rGrouped all {:,} ids into {:,} centers\n".format(i, c))
    return cids

def print_cids(cids, f):
    i = 1
    start = time.time()
    with open(f, 'w') as o:
        items = cids.items()
        t1, t2 = items.pop()
        o.write(">%s\n%s" %(t1, t2))
        for t1, t2 in items:
            o.write("\n>%s\n%s" %(t1, t2))
            i += 1
            if time.time() - start > 1 and args.verbose:
                sys.stdout.write("\rWrote {:,} entries to {}".format(i, f))
                sys.stdout.flush()
                start = time.time()
    if args.verbose:
        sys.stdout.write("\rWrote all {:,} entries to {}\n".format(i, f))

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
    m = []
    t = 0
    start = time.time()
    for center in centers:
        mhd = len(center)
        for cn in centers:
            dif = 0
            if cn != center:
                for a, b in zip(center, cn):
                    if a != b:
                        dif += 1
                    if dif >= mhd:
                        break
                mhd = dif
        m.append(mhd)
        t += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rFound Min Hamming Distance for {:,} centers".format(t))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rFound All Min Hamming Distance for {:,} centers\n".format(t))

    p = PdfPages(args.out.split('.')[0] + "_mhd.pdf")
    plot = plt.figure()
    plt.title("Minimum Hamming Distances Between Centers")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.xlabel("Minimum Hamming Distance")
    plt.ylabel("Centers")
    plt.hist(m, bins = range(0, len(centers[0])), align = 'left', color = 'blue', alpha = .6)
    plt.xlim(0, len(centers[0]))
    p.savefig(plot)
    p.close()

def jackpottogram(cids):
    n = [len(l) for l in cids.values()]
    total = sum(n)
    n = sorted(n, reverse = True)
    even = float(total) / len(cids)
    plot = plt.figure()

    slices = []
    labels = []
    rest = 0
    other = 0
    o = 0
    for rest in range(len(n)):
        per = n[rest] / float(total)
        if per > float(even)/2 or rest < 7:
            slices.append(per)
            labels.append('{:.3%} ({:,} reads)'.format(per, n[rest]))
        elif n[rest] == 1:
            break
        else:
            other += n[rest]
            o += 1
    per = other / float(total)
    slices.append(per)
    labels.append('Other {:,} Clusters:\n{:.3%} ({:,} reads)'.format(o, per, other))
    if rest + 1 != len(n):
        ones = (len(n)-rest)/float(total)
        slices.append(ones)
        labels.append("1's: {:.3%}".format(ones))

    color = cm.Pastel2(np.linspace(0.,1.,len(slices)))

    handles, text = plt.pie(slices, colors=color, startangle = 90)
    for handle in handles:
        handle.set_edgecolor('white')
        handle.set_linewidth(.05)
    plt.legend(handles, labels, title = "Largest Clusters", loc="upper right", prop={'size':8})
    plt.axis('equal')
    plt.title('Reads per Cluster by Percentage\nCut Barcodes <= {:,}'.format(args.cutreads))
    p = PdfPages(args.out.split('.')[0] + "_jkp.pdf")
    p.savefig(plot)
    p.close()

def rpc(cids):
    n = sorted([len(i) for i in cids.values()])

    #histogram binned by cluster size
    hist = plt.figure()
    plt.title("Reads per Cluster")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.xlabel("Cluster Size in Reads")
    plt.ylabel("Clusters")
    plt.hist(n, 155, color = 'blue', alpha = .6)

    #pie chart of number of different sized clusters
    slices, labels = pie_prep(n)
    pie = plt.figure()
    color = cm.Pastel2(np.linspace(0.,1.,len(slices)))
    handles, text = plt.pie(slices[::-1], colors=color, startangle = 90)
    for handle in handles:
        handle.set_edgecolor('white')
        handle.set_linewidth(.05)
    plt.legend(handles, labels[::-1], title = "Cluster Sizes", loc="upper right", 
        prop={'size':8})
    plt.axis('equal')
    plt.title('Number of Clusters at each Size')

    #cut outliers (1.5 * inter quartile range)
    if len(n) > 4:
        iqr = n[(3*len(n) + 3) / 4] - n[(len(n) + 1) / 4]
        lower = n[(len(n) + 1) / 4] - 1.5 * iqr
        upper = n[(3*len(n) + 3) / 4] + 1.5 * iqr
        f = [count for count in n if count > lower and count < upper]
    
        #histogram w/o outliers
        fhist = plt.figure()
        plt.title("Reads per Cluster w/o Outliers")
        plt.tick_params(axis = "both", labelsize = 8)
        plt.xlabel("Cluster Size in Reads between {:,} and {:,}".format(f[0], f[-1]))
        plt.ylabel("Clusters")
        plt.hist(f, 51, color = 'blue', alpha = .6)

        #cutoff pie chart
        slices, labels = pie_prep(f)
        fpie = plt.figure()
        color = cm.Pastel2(np.linspace(0.,1.,len(slices)))
        handles, text = plt.pie(slices[::-1], colors=color, startangle = 90)
        for handle in handles:
            handle.set_edgecolor('white')
            handle.set_linewidth(.05)
        plt.legend(handles, labels[::-1], title = "Cluster Sizes", loc="upper right", 
            prop={'size':8})
        plt.axis('equal')
        plt.title('Number of Clusters at each Size between {:,} and {:,}'.format(f[0], f[-1]))

    p = PdfPages(args.out.split('.')[0] + "_rpc.pdf")
    p.savefig(hist)
    if len(n) > 4:
        p.savefig(fhist)
    p.savefig(pie)
    if len(n) > 4:
        p.savefig(fpie)
    p.close()

def pie_prep(n):
    grp = 5
    tmp = [10000, 5000, 1000, 500, 250, 100, 50, 25, 10, 5]
    for x in range(n[-1]/27, n[-1]/5):
        for t in tmp:
            if x % t == 0:
                grp = x
                break
            elif grp % t == 0:
                break
    slices = []
    labels = []
    count = 0
    slc = 1
    for cl in n:
        if cl > grp * slc:
            slices.append(count)
            labels.append("{:,} - {:,} reads: {:,}".format(grp * (slc-1), grp * slc, count))
            count = 0
            slc += 1
        count += 1
    slices.append(count)
    labels.append("{:,} - {:,} reads: {:,}".format(grp * (slc-1), n[-1], count))
    return slices, labels

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
        of A, C, T, G will be considered.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("i", help = "fasta/q file containing sequences to be clustered")
    parser.add_argument("out", help = "output tag")
    parser.add_argument("-cr", "--cutreads", help = "cut out barcodes with <= reads than this value",
        type = float, default = 1)
    parser.add_argument("-cc", "--cutcollisions", help = "cut out clusters which probably include more\
        than one original barcode", action = "store_true", default = False)
    parser.add_argument("-g", "--graph", help = "creates graphs if specified", 
        action = "store_true", default = False)
    parser.add_argument("-v", "--verbose", help = "output progress information to terminal",
        action = "store_true", default = False)
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()