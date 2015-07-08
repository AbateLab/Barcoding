import os
import ast
import sys
import math
import time
import argparse
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#converts Q-score arguments to probabilities, and finishes
#argument parsing by actually calling the required functions.
def main():
    if args.eval:
        #qscores must be integers
        if args.basecutoff % 1 == 0:
            args.basecutoff = q2p(args.basecutoff)
        if args.basecomp % 1 == 0:
            args.basecomp = q2p(args.basecomp)
        if args.pcut % 1 == 0:
            args.pcut = q2p(args.pcut)
        evalq()
    if args.filter:
        filterq()
    if args.graph:
        seqEndN()
        baseComposition()
        scoreDistribution()
    if not(args.eval or args.filter or args.graph):
        print "no action requested"

#sums up cumulative probability incorrect of each read, places in folder/log
#replaces bases with prob wrong > "-b" with "N", puts reads in folder/nReplace
#for now also puts number of consecutive N's from the end in folder/log
def evalq():
    start = time.time()
    #check that the files and folders exist
    if os.path.exists(args.folder):
        print "warning: outputting into existing folder: " + args.folder
    else:
        os.makedirs(args.folder)
    if not os.path.isfile(args.eval):
        print "could not find input file: " + args.eval
        return

    with open(args.eval, 'r') as f, open(args.folder + "/log", "w") as log,\
    open(args.folder + "/nReplace", "w") as nRep, open(args.folder + "/report", "w") as rep:
        nCnt = [] #number of N's at each read position
        aCnt = [] #number of A's
        tCnt = [] #number of T's
        cCnt = [] #number of C's
        gCnt = [] #number of G's
        newN = [] #number of new N's
        t = 0 #total number of bases turned into N
        for lab, seq, exp in readfx(f):
            #actual summing loop, and replacing bases with N's and such
            e = 0 #expected number of bases incorrect in this read score
            m = 0 #max phred score
            p = 0 #min qualifying percent
            for i in range(len(exp)):
                #update those array things
                if seq[i].upper() == 'A':
                    while len(aCnt) <= i:
                        aCnt.append(0)
                    aCnt[i] += 1
                elif seq[i].upper() == 'T':
                    while len(tCnt) <= i:
                        tCnt.append(0)
                    tCnt[i] += 1
                elif seq[i].upper() == 'C':
                    while len(cCnt) <= i:
                        cCnt.append(0)
                    cCnt[i] += 1
                elif seq[i].upper() == 'G':
                    while len(gCnt) <= i:
                        gCnt.append(0)
                    gCnt[i] += 1
                elif seq[i].upper() == 'N':
                    while len(nCnt) <= i:
                        nCnt.append(0)
                    nCnt[i] += 1

                score = asc2p(exp[i])
                #find max phred score
                if score > m:
                    m = score
                #find percent > phred
                if score < args.pcut:
                    p += 1.0
                #find expected number of bases incorrect
                if score > args.basecutoff:
                    seq = list(seq)
                    og = seq[i].upper()
                    seq[i] = 'N'
                    seq = "".join(seq)
                    #if this read we're seeing is longer than any previous reads,
                    #and we need an N here, make the npos array longer
                    while len(newN) <= i:
                        newN.append(0)
                    if og != "N":
                        newN[i] += 1
                    if og == 'A':
                        aCnt[i] -= 1
                    elif og == "T":
                        tCnt[i] -= 1
                    elif og == "C":
                        cCnt[i] -= 1
                    elif og == "G":
                        gCnt[i] -= 1
                    e += args.basecomp
                else:
                    e += score
            p /= len(seq)
            i = -1 #seq[-1] is last character of sequence
            n = 0 #number of consecutive N's from the end of the sequence
            #start from the back, move forward. Really rough way of checking
            #whether Illumina's just filling in bases
            while n < len(seq) and seq[i] == 'N':
                n += 1
                i -= 1
            log.write(str(e) + " " + str(m) + " " + str(p) + " " +str(n) + "\n")
            nRep.write("@" + lab + "\n" + seq + "\n+\n" + exp + "\n")
            t += 1
            if t%10000 == 0 and args.verbose:
                sys.stdout.write("\rEval'd %i reads" %t)
                sys.stdout.flush()
        sys.stdout.write("\rEval'd all %i reads\n" %t)
        rep.write("Input file: " + args.eval + "\nBase Cutoff: " + str(args.basecutoff) +
            "\nEval runtime (min): " + str((time.time() - start) / 60) + "\nN Compensation: " +
            str(args.basecomp) + "\nnewN: " + str(newN) + "\naCnt: " + str(aCnt) + "\ntCnt: " +
            str(tCnt) + "\ncCnt: " + str(cCnt) + "\ngCnt: " + str(gCnt) + "\nnCnt: " + str(nCnt) +
            "\n---")

#splits folder/nReplace into folder/out and folder/out.cut. Out is what is good enough for you (-r)
#you can decide what's good enough by some of the graphs, and by trial and error
#Also removes end of reads that were decided to be "illumina zombie 'every read needs to be 
#the same length' filler reads"
def filterq():
    #check that files and folders exist
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/nReplace"):
        print args.folder + "/nReplace not found!\nDid you run the initial evaluation run?"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return

    start = time.time()
    scores = []
    with open(args.folder + "/nReplace", "r") as f, open(args.folder + "/log", "r") as log,\
        open(args.folder + "/" + args.filter, "w") as good,\
        open(args.folder + "/" + args.filter + ".cut", "w") as bad,\
        open(args.folder + "/report", "a") as rep:
        g = 0 #to report number of reads who made the cut
        b = 0 #to report number of reads who didn't
        cg = 0 #to report number of reads with N's cut from the back in the good pile
        cb = 0 #the above, for the bad
        t = 0
        for lab, seq, exp in readfx(f):
            l = log.next().split() #a line in folder/log
            n = int(l[-1])#the number of consecutive n's from the end
            #only do trimming if more n's in a row than ncutoff WILL THINK OF SOMETHING BETTER
            if n < args.ncutoff:
                n = -len(seq)
            if (args.fopt != 2 and float(l[args.fopt]) < args.readcutoff) or (args.fopt == 2 and 
                    float(l[args.fopt]) > args.readcutoff) and n != len(seq):
                scores += [float(l[args.fopt])]
                good.write("@" + lab + "\n" + seq[:-n] + "\n+\n" + exp[:-n]+"\n")
                if n >= args.ncutoff:
                    cg += 1
                g += 1
            else:
                bad.write("@" + lab + "\n" + seq[:-n] + "\n+\n" + exp[:-n]+"\n")
                if n >= args.ncutoff:
                    cb += 1
                b += 1
            t += 1
            if t%50000 == 0 and args.verbose:
                sys.stdout.write("\rFltr'd %i reads" %t)
                sys.stdout.flush()
        sys.stdout.write("\rFltr'd all %i reads\n" %t)
        fopts = ["Expected Incorrect Bases", "Max Allowable Phred", "Percent Qualifying Bases"]
        rep.write("\n\nFilter to '" + str(args.filter) + "'\nFilter option: " + fopts[args.fopt] +
            "\nFilter runtime (min): " + str((time.time() - start) / 60) +
            "\nRead Error Cutoff: " + str(args.readcutoff) +
            "\nSequential Ending N Cutoff: " + str(args.ncutoff) +
            "\nPercent Reads Retained: " + str(float(g) / (g + b)) +
            "\nNumber of Remaining Reads: " + str(g) + "\n\tNumber of Trimmed Reads: " + str(cg) +
            "\nNumber of Cut Reads: " + str(b) + "\n\tNumber of Trimmed Reads: " + str(cb))
    if scores:
        sHist(scores, ["Expected Incorrect Bases per Read", "Maximum Single Base Error", 
            "Percent of Read Better Than pcut"][args.fopt], False, args.filter)
    else:
        print "Cutoff too stringent - no reads made the cut"

#run through log to gather info, spit out a graph
def seqEndN():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return
    n = []
    m = 0
    with open(args.folder + "/log", "r") as log:
        for line in log:
            l = line.split()
            if int(l[-1]) > m:
                m = int(l[-1])
            n.append(int(l[-1]))
    p = PdfPages(args.folder + "/seqN.pdf")
    plot = plt.figure()
    plt.title("cumulative Reads Ending in x Sequential N's")
    plt.tick_params(axis = "both", labelsize = 8)
    plt.ylabel("Number of Sequential N's from end")
    plt.xlabel("cumulative Number of Reads")
    plt.hist(n, 257, color = 'blue', histtype = 'step', orientation = 'horizontal', 
        cumulative=True, fill = True, alpha = .3)
    plt.ylim([-.1, m + .1])
    plt.gca().invert_yaxis()
    p.savefig(plot)
    p.close()    

#make graph of # of actg's at each position
def baseComposition():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/report"):
        print args.folder + "/report not found!\nDid you run the initial evaluation run?"
        return

    a = []
    t = []
    c = []
    g = []
    n = []
    nn = []
    with open(args.folder + "/report", "r") as rep:
        for line in rep:
            if line == "---\n":
                break
            lab, dat = line.split(": ")
            if lab == "newN":
                nn = ast.literal_eval(dat[:-1])
            elif lab == "aCnt":
                a = ast.literal_eval(dat[:-1])
            elif lab == "tCnt":
                t = ast.literal_eval(dat[:-1])
            elif lab == "cCnt":
                c = ast.literal_eval(dat[:-1])
            elif lab == "gCnt":
                g = ast.literal_eval(dat[:-1])
            elif lab == "nCnt":
                n = ast.literal_eval(dat[:-1])
    #make everything equal length
    m = max(len(i) for i in [a, t, c, g, n, nn])
    for l in [a, t, c, g, n, nn]:
        if len(l) < m:
            l.extend([0]*(m - len(l)))
    tot = [sum(vals) for vals in zip(a, t, c, g, n, nn)]
    x = list(range(m))

    p = PdfPages(args.folder + "/baseContent.pdf")
    plot = plt.figure()
    ax = plt.subplot(111)
    ax.set_title("Bases at each Read Position")
    ax.tick_params(axis = "both", labelsize = 8)
    ax.set_ylabel("Bases")
    ax.set_xlabel("Position")
    ax.plot(x, a, alpha = .9, linewidth = 1)
    ax.plot(x, t, alpha = .85, linewidth = 1)
    ax.plot(x, c, alpha = .8, linewidth = 1)
    ax.plot(x, g, alpha = .75, linewidth = 1)
    ax.plot(x, n, alpha = .7, linewidth = 1)
    ax.plot(x, nn, alpha = .65, linewidth = 1)
    ax.plot(x, tot, linewidth = 1)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(['A', "T", "C", "G", "Existing N", "Added N", "Total"], loc = "center left",
        bbox_to_anchor=(1, .5), fancybox = True, ncol = 1)
    p.savefig(plot)
    p.close()

def scoreDistribution():
    if not os.path.exists(args.folder):
        print folder + " not found!"
        return
    elif not os.path.isfile(args.folder + "/log"):
        print args.folder + "/log not found!\nDid you run the initial evaluation run?"
        return
    exp = []
    mxp = []
    pcp = []
    with open(args.folder + "/log", "r") as log:
        for line in log:
            e, m, p, n = line.split()
            exp += [float(e)]
            mxp += [float(m)]
            pcp += [1-float(p)]
    sHist(exp, "Expected Incorrect Bases per Read", True, 'expHist')
    sHist(mxp, "Maximum Single Base Error", True, 'mxpHist')
    sHist(pcp, "Percent of Read Worse Than pcut", True, 'pcpHist')

def sHist(scores, heuristic, cum, fname):
    p = PdfPages(args.folder + "/" + fname + ".pdf")
    plot = plt.figure()
    scores.sort()
    plt.hist(scores, bins = 251, cumulative = cum, histtype = 'step', fill = True, alpha = .3)
    plt.title(heuristic + "%s Histogram" %" Cumulative" if cum else "")
    plt.xlabel(heuristic)
    plt.ylabel("Number of Reads")
    p.savefig(plot)
    p.close()

def asc2p(asc):
    return q2p(asc2q(asc))

def asc2q(asc):
    return ord(asc) - 33

def p2q(p):
    return abs(10 * math.log(p, 10))

def q2p(q):
    return 10**(-float(q) / 10)

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

#Argument Parsing, help text, and calling of the main function if enough arguments present
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="qScoreQC evaluates and filters phred fastq files.\
        One must run 'eval' once before filtering and graphing, but may run other options\
        limitlessly afterwards. One can run everything together, qScoreQC will simply call eval first.\
        A basic run looks like: python qScoreQC.py myfolder -e -i input_file -g -f -o output_tag.\
        This evals, graphs, and filters, all with default values. After eval, qScoreQC provides five\
        different graphs - number of each base at each position, N tail lengths cumulative histogram,\
        Expected Errors per Read histogram, Min Q-score per Read histogram, Percent Read above Q-Score\
        Histogram - to assist the selection of filtering parameters. You may filter by any of three\
        heuristics - Errors per Read (summation of probabilities associated with Q-score), minimum\
        allowable Q-score per read, and minimum percent of read above selected Q-score", 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("folder", help = "folder outputting to, or working in")
    parser.add_argument("-e", "--eval", help = "run Eval function on specified fastq file", default = "")
    parser.add_argument("-pc", "--pcut", help = "cutoff for quality of base used in the percent qualifying\
        filter. Bases with Q-Score better than pcut count towards 'percent qualifying'.\
        Interprets as phred Q-score if integral, and as probability incorrect otherwise.", 
        default = 30, type = float)
    parser.add_argument("-b", "--basecutoff", help = "cutoff for base conversion to N. Bases with Q-score\
        worse than basecutoff are replaced with N's. Interprets as phred Q-score if integral, and as\
        probability incorrect otherwise.", default = 20, type = float)
    parser.add_argument("-bc", "--basecomp", help = "reported Q-score of a base replaced by N,\
        interprets as phred Q-score if integral, and as probability incorrect otherwise.", 
        default = 22, type = float)

    parser.add_argument("-f", "--filter", help = "run Filter function, outputting using specified tag",
        default = "")
    parser.add_argument("-fo", "--fopt", help = "select filtering heuristic with 0: Expected Number of\
        Incorrect Bases per Read, 1: Maximum Allowable Q-score, 2: Percent of Read over pcut", 
        type=int, choices = [0, 1, 2], default = 0)
    parser.add_argument("-n", "--ncutoff", help = "Reads which end in over ncutoff sequential N's will\
        have all sequential N's from end trimmed", default = 1000, type = int)
    parser.add_argument("-r", "--readcutoff", help = "fo:0 - Reads with Exp score > readcutoff will be\
        cut (0 to len(read)), fo:1 - Reads whose worst base is worse than readcutoff will be cut\
        interprets as Q-score if integral, and probability incorrect otherwise, fo:2 - Reads with less\
        than readcutoff percent of bases better than pcut will be cut.", default = .001, type = float)
    
    parser.add_argument("-g", "--graph", help = "creates graphs if specified", 
        action = "store_true", default = False)
    parser.add_argument("-v", "--verbose", help = "output progress information to terminal",
        action = "store_true", default = False)
    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        main()