import sys
import time
import math
import random
import argparse
import matplotlib
matplotlib.use('Agg')
import operator as op
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
    ccounts = initialize(args.len, args.num)
    if args.graph:
        minDistHist(ccounts.keys())
    t = sum(ccounts.values())
    u = len(ccounts.keys())
    l = len(ccounts.keys()[0])
    j = 0
    for i in range(args.cyc):
        ccounts, t, u, j = stats_cyc(ccounts, args.err, args.cpr, i+1, t, u, j, l)
    ccounts = sample(ccounts, args.sample, t)
    cids = ccounts2cids(ccounts)
    cids2fasta(cids, args.out+".fasta")

def initialize(l, n):
    ccounts = {}
    i = 0
    start = time.time()
    for _ in range(n):
        ccounts[''.join(random.SystemRandom().choice('ACTG') for _ in range(l))] = 1
        i += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rInitialized %i Barcodes" %i)
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rInitialized all %i Barcodes\n" %i)
    if len(ccounts) < n:
        print "Barcode Collision Detected"
    with open(args.out + "_ocn.txt", 'w') as f:
        keys = ccounts.keys()
        f.write(keys.pop())
        for barcode in keys:
            f.write("\n" + barcode)
    return ccounts

def stats_cyc(ccounts, er, cr, cyc, t, u, j, l):
    c = ccounts.copy()
    start = time.time()
    for center, count in c.items():
        if count > args.lnlaw:
            rem, count = math.modf(count * cr)
            count = int(count)
            if random.random() < rem:
                count += 1
            n = (1 - er)**l
            rem, num = math.modf(count * n)
            num = int(num)
            if random.random() < rem:
                num += 1
            c[center] += num
            t += num
            count -= num
            weight = 1 - n
            num_err = 1
            while count > 0:
                n = ncr(l, num_err) * (er**num_err) * ((1 - er)**(l - num_err))
                rem, num = math.modf(count * n / weight)
                num = int(num)
                weight -= n
                if random.random() < rem:
                    num += 1
                for _ in range(num):
                    new = list(center)
                    for i in random.sample(range(l), num_err):
                        chars = ['A', 'T', 'C', 'G']
                        chars.remove(new[i].upper())
                        new[i] = random.choice(chars)
                    new = ''.join(new)
                    if new in c:
                        c[new] += 1
                    else:
                        c[new] = 1
                        u += 1
                    if num_err >= 2:
                        j += 1
                    t += 1
                count -= num
                num_err += 1 if num_err < l else 0
        else:
            for _ in range(count):
                if random.random() < cr:
                    new = []
                    d = 0
                    for char in center:
                        if random.random() < args.err:
                            chars = ['A', 'T', 'C', 'G']
                            chars.remove(char)
                            new.append(random.choice(chars))
                            d += 1
                        else:
                            new.append(char)
                    new = ''.join(new)
                    if new in c:
                        c[new] += 1
                    else:
                        c[new] = 1
                        u += 1
                    if d >= 2:
                        j += 1
                    t += 1
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rCycle: %i\tJumps: %i\tUnique: %i\tTotal: %i" %(cyc, j, u, t))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rCycle: %i\tJumps: %i\tUnique: %i\tTotal: %i\n" %(cyc, j, u, t))
    return c, t, u, j

def sample(ccounts, num, tot):
    index = 0
    c = {}
    if num >= tot:
        return ccounts
    samp = random.sample(xrange(tot), num)
    samp.sort(reverse = True)
    start = time.time()
    items = ccounts.items()
    tmp = samp.pop()
    for center, count in items:
        if not samp:
            break
        index += count
        i = 0
        while tmp < index:
            i += 1
            if samp:
                tmp = samp.pop()
            else:
                break
        if i > 0:
            c[center] = i
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rSampled %i reads" %(num - len(samp)))
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rSampled %i reads\n" %(num - len(samp)))
    return c


def ccounts2cids(ccounts):
    cids = {}
    i = 0
    start = time.time()
    for center, count in ccounts.items():
        ids = []
        for _ in range(count):
            ids.append(i)
            i += 1
        cids[center] = ids
        if time.time() - start > 1 and args.verbose:
            sys.stdout.write("\rId'd %i sequences" %i)
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rId'd all %i sequences\n" %i)
    return cids

def cids2fasta(cids, f):
    start = time.time()
    with open(f, 'w') as f:
        items = cids.items()
        center, i = items[0][0], items[0][1].pop()
        f.write(">%i\n%s" %(i, center))
        n = 1
        for center, ids in items:
            for i in ids:
                f.write("\n>%i\n%s" %(i, center))
                n += 1
            if time.time() - start > 1 and args.verbose:
                sys.stdout.write("\rWrote %i sequences" %n)
                sys.stdout.flush()
                start = time.time()
    if args.verbose:
        sys.stdout.write("\rWrote all %i sequences\n" %n)

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
            sys.stdout.write("\rFound Min Hamming Distance for %i Barcodes" %t)
            sys.stdout.flush()
            start = time.time()
    if args.verbose:
        sys.stdout.write("\rFound All Min Hamming Distance for %i Barcodes\n" %t)

    p = PdfPages(args.out.split('.')[0] + "_mhd.pdf")
    plot = plt.figure()
    plt.title("Minimum Hamming Distances Between Barcodes")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.xlabel("Minimum Hamming Distance")
    plt.ylabel("Barcodes")
    plt.hist(m, bins = range(0, len(centers[0])), align = 'left', color = 'blue', alpha = .6)
    plt.xlim(0, len(centers[0]))
    p.savefig(plot)
    p.close()

def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Creates Fasta File of simulated PCR run on barcodes", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("out", help = "output fasta file will be named [out].fasta, output minimum\
        hamming distance histogram will be named [out]_mdh.pdf, and original barcodes will be placed in\
        [out].txt")
    parser.add_argument("-l", "--len", help = "length of barcodes", default = 15, type = int)
    parser.add_argument("-n", "--num", help = "number of initial barcodes", default = 10000, type = int)
    parser.add_argument("-e", "--err", help = "probability a base is copied incorrectly",
        type = float, default = .0001)
    parser.add_argument("-r", "--cpr", help = "probability a sequence is copied", type = float, 
        default = .8)
    parser.add_argument("-c", "--cyc", help = "number of simulated pcr cycles", type = int, 
        default = 25)
    parser.add_argument("-lln", "--lnlaw", help = "number of id's per barcode above which simulation\
        starts using statistics to speed up simulation. Will automatically partition away expected\
        number of sequences to be copied, copied without error, with one error, so on", type = int,
        default = 1000)
    parser.add_argument("-s", "--sample", help = "number of sequences to be kept for output", type = int,
        default = 10000000)
    parser.add_argument("-g", "--graph", help = "creates graphs if specified", 
        action = "store_true", default = False)
    parser.add_argument("-v", "--verbose", help = "output progress information to terminal",
        action = "store_true", default = False)
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()