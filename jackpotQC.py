import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import numpy as np
import random

def main():
    counts = sorted(make_set(args.i).values(), reverse = True)
    #report(counts)
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
    print "max: " + str(counts[0]) + "\nmin: " + str(counts[-1])
    total = 0
    for x in range(len(counts)):
        total += counts[x]
    for x in range(0,5):
        per = 100*(counts[x] / float(total))
        perstr = '%.2f' % per
        print str(counts[x]) + ' reads: ' + perstr + '%'

def jackpottogram(counts):
    fig1 = plt.figure() #cumulative histogram
    counts.reverse() #make into ascending
    cumulvals = []
    cumulvals.append(counts[0])
    for x in range(1,len(counts)):
        cumulvals.append(counts[x] + cumulvals[x-1])
    x_range = range(len(counts))
    n = len(counts)
    total = float(cumulvals[n-1])
    slope = total / n
    expected = []
    for i in x_range:
        expected.append(i*slope)
    x, y = lowres(x_range, cumulvals)
    plt.plot(x, y)
    plt.fill_between(x, y, 0, color='blue', alpha=0.5)
    x, y = lowres(x_range, expected)
    plt.plot(x, y, 'r--') #expected distr 
    plt.plot(x, x, 'r--') #slope of 1
    plt.title('Cumulative Histogram of Barcode Reads')
    plt.xlabel('Barcodes')
    plt.ylabel('Reads')
    plt.tick_params(axis='both', which='major', labelsize='8')
    plt.annotate('Slope for even distribution(' + str(int(slope)) + ')', xy=(x_range[n-100000], 
        expected[n-100000]), xytext=(x_range[n-1]-500000, expected[n-1]-100000), 
        arrowprops=dict(facecolor='black', shrink=0.05),
        horizontalalignment='center', fontsize='8') #wont look nice for every graph u plot... v arbritrary
    plt.annotate('Slope = 1', xy=(x_range[n-300000], x_range[n-300000]), 
        xytext=(x_range[n-1]-300000, x_range[n-1]+100000), arrowprops=dict(facecolor='black', shrink=0.05),
        horizontalalignment='left', fontsize='8')
    p = PdfPages(args.o + "_jckptogram.pdf")
    p.savefig(fig1)
    fig2 = plt.figure() #pie chart
    counts.reverse() #make descending
    pervals = []
    labels = []
    rest = 0
    for x in range(100): #first 100 gets own slice
        per = 100*(counts[x]/float(total))
        if per < 1.0 and x > 10:
            rest = x
            break
        else:
            pervals.append(per)
            perstr = '%.3f' % per
            labels.append(perstr + '%')
    other = 0
    for x in range(rest, len(counts)): #rest are grouped into own slice
        other += counts[x]
    other = 100* (other/float(total))
    pervals.append(other)
    otherstr = 'Other: %.3f' % other
    labels.append(otherstr + '%')
    ncolors = rest + 1
    cs = [] #for a rainbow gradient of colors
    color = iter(cm.rainbow(np.linspace(0,1,ncolors)))
    for i in range(ncolors-1):
        c = next(color)
        cs.append(c)
    random.shuffle(cs)
    cs.append('grey') #other slice = always grey
    handles, text = plt.pie(pervals, colors=cs, startangle = 90)
    for handle in handles:
        handle.set_edgecolor('white')
        handle.set_linewidth(.05)
    plt.legend(handles, labels, title = "Largest Barcodes", loc="upper right", prop={'size':8})
    plt.axis('equal')
    plt.title('Reads per Barcode by Percentage')
    p.savefig(fig2)
    p.close()

def lowres(x, y):
    t = 1
    if len(x) != len(y):
        return "lists must be same length"
    while len(x) / t > 1000:
        t += 1
    lx = []
    ly = []
    for i in range(0, len(x), t):
        lx.append(x[i])
        ly.append(y[i])
    if i + 1 < len(x):
        lx.append(x[-1])
        ly.append(y[-1])
    return lx, ly

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