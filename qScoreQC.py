import sys
import os
import math
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def main():
    if args.eval:
        if args.input:
            if args.basecutoff % 1 == 0:
                args.basecutoff = q2p(args.basecutoff)
            if args.basecomp % 1 == 0:
                args.basecomp = q2p(args.basecomp)
            evalq()
        else:
            print "eval requires -i to run. -b defaults to 20, -bc defaults to 25"
            return
    if args.filter:
        if args.out:
            filterq()
        else:
            print "filter requires -o to run. -n defaults to 1000, -r defaults to .001"
    if args.fracError:
        fracErrorGraph()
    if args.cummToler:
        cummToleranceGraph()
    if args.nPosDistr:
        nPosDistribution()
    if args.seqN:
        seqEndN()
    if not(args.eval or args.filter or args.fracError or args.cummToler or args.nPosDistr or args.seqN):
        print "no action requested"

def evalq():
    if os.path.exists(args.folder):
        print "warning: outputting into existing folder: " + args.folder
    else:
        os.makedirs(args.folder)
    if not os.path.isfile(args.input):
        print "could not find input file: " + args.input
        return
    with open(args.input, 'r') as f, open(args.folder + "/log", "w") as log,\
    open(args.folder + "/nReplace", "w") as nRep:
        f.next()
        npos = [0] * len(f.next().rstrip('\n'))
        f.seek(0)
        for line in f:
            lab = line
            seq = f.next()
            pls = f.next()
            exp = f.next()
            if len(seq.rstrip('\n')) != len(exp.rstrip('\n')):
                print "sequence and qscore different lengths: " + lab
                return
            e = 0
            n = 0
            i = -2            
            for i in range(len(exp)-1):
                score = asc2p(exp[i])
                if score > args.basecutoff or seq[i] == 'N':
                    seq = list(seq)
                    seq[i] = 'N'
                    seq = "".join(seq)
                    npos[i] += 1.0
                    e += args.basecomp
                else:
                    e += score
            while seq[i] == 'N' and n < len(seq):
                n += 1
                i -= 1
            log.write(str(e) + " " + str(n) + "\n")
            nRep.write(lab)
            nRep.write(seq)
            nRep.write(pls)
            nRep.write(exp)
        log.write(str(npos))

def filterq():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/nReplace"):
        print args.folder + "/nReplace not found!\nDid you run the initial evaluation run?"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return
    with open(args.folder + "/nReplace", "r") as f, open(args.folder + "/log", "r") as log,\
        open(args.folder + "/" + args.out, "w") as good,\
        open(args.folder + "/" + args.out + ".cut", "w") as bad:
        for line in f:
            l = log.next().split()
            n = int(l[1])
            if n < args.ncutoff:
                n = 0
            if float(l[0]) < args.readcutoff:
                good.write(line)
                good.write(f.next()[:-n-1]+"\n")
                good.write(f.next())
                good.write(f.next()[:-n-1]+"\n")
            else:
                bad.write(line)
                bad.write(f.next()[:-n-1]+"\n")
                bad.write(f.next())
                bad.write(f.next()[:-n-1]+"\n")                

def fracErrorGraph():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return
    total = 0
    noerr = 0
    with open(args.folder + "/log", "r") as log:
        for line in log:
            l = line.split()
            if l[0][0] != '[':
                if float(l[0]) < args.fracErrorCutoff:
                    noerr += 1
                total += 1
    p = PdfPages(args.folder + "/frcError.pdf")
    t = list(xrange(1,total+1))
    s = [min(1, noerr / float(n)) for n in t]
    plot = plt.figure()
    plt.title("%Reads with Expected Number of Incorrect Bases < " + 
        str(args.fracErrorCutoff) + "\n" +str(noerr) + " Reads < " + str(args.fracErrorCutoff))
    plt.ylim([-.01, 1.01])
    plt.ylabel("Fraction of Reads Under Threshold")
    plt.xlabel("Number of Reads Used")
    plt.plot(t, s)
    p.savefig(plot)
    p.close()

def cummToleranceGraph():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return
    exp = []
    m = 0
    with open(args.folder + "/log", "r") as log:
        for line in log:
            l = line.split()
            if l[0][0] != '[':
                if float(l[0]) > m:
                    m = float(l[0])
                exp.append(float(l[0]))
    p = PdfPages(args.folder + "/expHist.pdf")
    plot = plt.figure()
    plt.title("Cummulative Reads by Quality")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.ylabel("Expected Number of Errors per Read")
    plt.xlabel("Cummulative Number of Reads")
    plt.hist(exp, 257, color = 'blue', histtype = 'step', orientation = 'horizontal', 
        cumulative=True, alpha = .6)
    plt.ylim([-.01, m + .01])
    plt.gca().invert_yaxis()
    p.savefig(plot)
    p.close()

def nPosDistribution():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return
    npos = None
    t = 0
    with open(args.folder + "/log", "r") as log:
        for line in log:
            if line[0] == '[':
                npos = [float(string) for string in line[1:-1].split(", ")]
            else:
                t += 1
    if npos == None:
        print "Could not find n position distribution information in log. Rerun eval?"
        return
    else:
        npos = [num / t for num in npos]
    p = PdfPages(args.folder + "/nDistr.pdf")
    plot = plt.figure()
    plt.tick_params(axis = 'both', labelsize = 8)
    plt.ylabel("Fraction of Reads with N")
    plt.xlabel("Position on Read")
    plt.title("%Reads with N's at each Position")
    plt.plot(list(range(len(npos))), npos, color='blue', linestyle='dashed', 
        marker='o',markerfacecolor='cyan', markersize=6, alpha = .6)
    p.savefig(plot)
    p.close()

def seqEndN():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return    
    exp = []
    m = 0
    with open(args.folder + "/log", "r") as log:
        for line in log:
            l = line.split()
            if l[0][0] != '[':
                if int(l[1]) > m:
                    m = int(l[1])
                exp.append(int(l[1]))
    p = PdfPages(args.folder + "/seqN.pdf")
    plot = plt.figure()
    plt.title("Cummulative Reads Ending in x Sequential N's")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.ylabel("Number of Sequential N's from end")
    plt.xlabel("Cummulative Number of Reads")
    plt.hist(exp, 257, color = 'blue', histtype = 'step', orientation = 'horizontal', 
        cumulative=True, alpha = .6)
    plt.ylim([-.1, m + .1])
    plt.gca().invert_yaxis()
    p.savefig(plot)
    p.close()    

def asc2p(asc):
    return q2p(asc2q(asc))

def asc2q(asc):
    return ord(asc) - 33

def p2q(p):
    return -10 * math.log(p, 10)

def q2p(q):
    return 10**(-float(q) / 10)

def roughQ2P(q):
    p = [1.0, 0.7943282347242815, 0.6309573444801932, 0.5011872336272722,
    0.3981071705534972, 0.31622776601683794, 0.251188643150958,
    0.19952623149688797, 0.15848931924611134, 0.12589254117941673,
    0.1, 0.07943282347242814, 0.06309573444801933, 0.05011872336272722,
    0.039810717055349734, 0.03162277660168379, 0.025118864315095794,
    0.0199526231496888, 0.015848931924611134, 0.012589254117941675, 0.01,
    0.007943282347242814, 0.00630957344480193, 0.005011872336272725,
    0.003981071705534973, 0.0031622776601683794, 0.0025118864315095794,
    0.001995262314968879, 0.001584893192461114, 0.0012589254117941675,
    0.001, 0.0007943282347242813, 0.000630957344480193, 0.0005011872336272725,
    0.00039810717055349735, 0.00031622776601683794, 0.00025118864315095795,
    0.00019952623149688788, 0.00015848931924611142, 0.00012589254117941674,
    0.0001, 7.943282347242822e-05, 6.309573444801929e-05, 5.011872336272725e-05,
    3.9810717055349695e-05, 3.1622776601683795e-05, 2.5118864315095822e-05,
    1.9952623149688786e-05, 1.584893192461114e-05, 1.2589254117941661e-05,
    1e-05, 7.943282347242822e-06, 6.30957344480193e-06, 5.011872336272725e-06,
    3.981071705534969e-06, 3.162277660168379e-06, 2.5118864315095823e-06,
    1.9952623149688787e-06, 1.584893192461114e-06, 1.2589254117941661e-06,
    1e-06, 7.943282347242822e-07, 6.30957344480193e-07, 5.011872336272725e-07,
    3.981071705534969e-07, 3.162277660168379e-07, 2.5118864315095823e-07,
    1.9952623149688787e-07, 1.584893192461114e-07, 1.2589254117941662e-07,
    1e-07, 7.943282347242822e-08, 6.30957344480193e-08, 5.011872336272725e-08,
    3.981071705534969e-08, 3.162277660168379e-08, 2.511886431509582e-08,
    1.9952623149688786e-08, 1.5848931924611143e-08, 1.2589254117941661e-08,
    1e-08, 7.943282347242822e-09, 6.309573444801943e-09, 5.011872336272715e-09,
    3.981071705534969e-09, 3.1622776601683795e-09, 2.511886431509582e-09,
    1.9952623149688828e-09, 1.584893192461111e-09, 1.2589254117941663e-09,
    1e-09, 7.943282347242822e-10, 6.309573444801942e-10, 5.011872336272714e-10,
    3.9810717055349694e-10, 3.1622776601683795e-10, 2.511886431509582e-10,
    1.9952623149688828e-10, 1.584893192461111e-10, 1.2589254117941662e-10]
    return p[q]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="qScoreQC evaluates and filters phred fastq files using\
        q-scores. One must run 'eval' once before running other options, but may run other options\
        limitlessly afterwards. One can run everything together, qScoreQC will simply call eval first.\
        When running 'eval', one must also  provide '-i', '-b', '-bc'. When running 'filter', one must\
        also provide '-o', 'n', 'r'.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder", help = "folder outputting to, or working in")
    parser.add_argument("-e", "--eval", help = "run Eval function", action = "store_true")
    parser.add_argument("-i", "--input", help = "fastq file to be evaluated")
    parser.add_argument("-b", "--basecutoff", help = "cutoff for base conversion to N, \
        interprets as phred q-score if integral, and as probability incorrect otherwise.", 
        default = 20, type = float)
    parser.add_argument("-bc", "--basecomp", help = "reported error level of a base replaced by N,\
        interprets as phred q-score if integral, and as probability incorrect otherwise.", 
        default = 25, type = float)
    parser.add_argument("-f", "--filter", help = "run Filter function", action = "store_true")
    parser.add_argument("-o", "--out", help = "fastq file of filtered reads")
    parser.add_argument("-n", "--ncutoff", help = "if read ends in over x number of N's, all sequential\
        N's at end of read will be disarded. WARNING: This may create a file with unequal read lengths!",
        default = 1000, type = int)
    parser.add_argument("-r", "--readcutoff", help = "expected number of bases incorrect per read\
        above which read will be filtered out.", default = .001, type = float)
    parser.add_argument("-fe", "--fracError", help = "graph fraction reads with error", 
        action = "store_true")
    parser.add_argument("-fec", "--fracErrorCutoff", help = "specify error level for fracError graph",
        default = 1, type = float)
    parser.add_argument("-ct", "--cummToler", help = "graph number of reads < some error rate", 
        action = "store_true")
    parser.add_argument("-nd", "--nPosDistr", help = "graph fraction of N's per position on read", 
        action = "store_true")
    parser.add_argument("-sn", "--seqN", help = "graph number of reads ending in x sequential N's",
        action = "store_true")
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()
