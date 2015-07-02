import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
    counts = None
    counts = make_set(args.i)
    counts = sorted(counts.values(), reverse = True)
    report(counts)
    jackpottogram(counts)

#reads fasta file to get set of barcodes to cluster
def make_set(f):
    counts = {}
    with open(f, 'r') as f:
        for lab, seq, exp in readfx(f):
            if seq in counts:
                counts[seq] += 1
            else:
                counts[seq] = 1
    return counts

def report(counts):
    n = 0
    print "max: " + str(counts[0]) + "\nmin: " + str(counts[-1])
    c = 10**(len(str(counts[0]))-1)
    for x in counts:
        if x < c:
            print str(n) + "\tbarcodes where " + str(c) + "\t<= reads/barcode < " + str(c*10)
            c /= 10
            n = 0
        n += 1
    print str(n) + "\tbarcodes where " + str(c) + "\t<= reads/barcode < " + str(c*10)

def jackpottogram(counts):
    fig = plt.figure()
    #TODO: plot jackpottogram
    p = PdfPages(args.o + "_jckptogram.pdf")
    p.savefig(fig)
    p.close()

def lowres(counts):
    low = []
    res = counts[0] / 100.
    clear = True
    tmp = 0
    for x in counts:
        if clear:
            if x < res:
                clear = False
                tmp += x
            else:
                low.append(x)
        else:
            if tmp >= res:
                low.append(tmp)
            else:
                tmp += x
    return low

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
    parser = argparse.ArgumentParser(description="summarizes jackpotting level of pcr results")
    parser.add_argument("i", help = "fasta/q formatted input")
    parser.add_argument("o", help = "tag for output files")
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()