import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
    counts = None
    counts = make_set(args.i)
    countslist = sorted(counts.values(), reverse = True)
    report(countslist)
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
    fig = plt.figure()
    #TODO: plot jackpottogram
    histvals = []
    for val in sorted(counts, key=counts.get): #sorts dict key by ascending frequency
        histvals.append(counts[val])
    cumulvals = [] #now make it cumulative
    cumulvals.append(histvals[0]) #first val is just the smallest val
    for x in range(1,len(histvals)): #x represents indexes from 1 to n
        cumulvals.append(histvals[x] + cumulvals[x-1]) #successive vals are cumulative    
    x_range = range(len(counts))
    n = len(counts)
    slope = float(cumulvals[n-1]) / n
    expected = []
    for i in x_range:
        expected.append(i*slope)
    plt.plot(x_range, cumulvals)
    plt.fill_between(x_range, cumulvals, 0, color='blue', alpha=0.5)
    plt.plot(x_range, expected, 'r--') #expected distr
    plt.plot(x_range, x_range, 'r--') #slope of 1
    plt.title('Cumulative Histogram of Barcode Reads')
    plt.xlabel('Barcodes')
    plt.ylabel('Reads')
    plt.tick_params(axis='both', which='major', labelsize='8')
    plt.annotate('Slope for even distribution(' + str(int(slope)) + ')', xy=(x_range[n-100000], expected[n-100000]), 
        xytext=(x_range[n-1]-500000, expected[n-1]-100000), arrowprops=dict(facecolor='black', shrink=0.05),
        horizontalalignment='center', fontsize='8') #wont look nice for every graph u plot... v arbritrary
    plt.annotate('Slope = 1', xy=(x_range[n-300000], x_range[n-300000]), 
        xytext=(x_range[n-1]-300000, x_range[n-1]+100000), arrowprops=dict(facecolor='black', shrink=0.05),
        horizontalalignment='left', fontsize='8')
    plt.show()
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