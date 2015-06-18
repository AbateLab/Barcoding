import sys
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def hist():
    l = []
    with open(args.i, 'r') as f:
        for line in f:
            l.append(len(f.next().rstrip("\n")))
            f.next()
            f.next()
    p = PdfPages(args.o + "RLHist.pdf")
    plot = plt.figure()
    plt.title("Read Length")
    plt.ylabel("Number of Reads")
    plt.xlabel("Read Length")
    plt.hist(l, 257, color = 'blue', alpha = .6)
    p.savefig(plot)
    p.close()    

parser = argparse.ArgumentParser(description = "plots a histogram of read lenghts from a fastq file into a pdf")
parser.add_argument("i", help = "input fastq file")
parser.add_argument("o", help = "output pdf destination")
if len(sys.argv) < 2:
    parser.print_help()
else:
    args = parser.parse_args()
    hist()