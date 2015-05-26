import sys
import os
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


### Global Default Values ###
inf = None
out = None
folder = None
rcut = .001
rGraph = True
bcut = .99
bcutcomp = .05
bGraph = True

def main():
    args = sys.argv
    global inf
    global out
    global folder
    global rcut
    global rGraph
    global bcut
    global bcutcomp
    global bGraph
    if args[1] == '--e':
        inf = args[2]
        args = args[3:]
        while args:
            opt = args.pop(0)
            if len(args) == 0:
                print "missing argument after " + opt
                help()
                return
            if opt == '-f':
                folder = args.pop(0)
            elif opt == '-bq':
                bcut = float(args.pop(0))
            elif opt == '-bc':
                bcutcomp = float(args.pop(0))
            elif opt == '-nb':
                bGraph = False
            elif opt == '-nr':
                rGraph = False
            else:
                print "nonsensical argument: " + opt
                help()
                return
        if folder == None:
            folder = inf + '_' + str(bcut).split('.')[1] + '_' + str(bcutcomp).split('.')[1]
        evalq()
    elif args[1] == '--f':
        folder = args[2]
        args = args[3:]
        while args:
            opt = args.pop(0)
            if opt == '-o':
                out = args.pop(0)
            elif opt == '-rc':
                rcut = float(args.pop(0))
            else:
                print "nonsensical argument: " + opt
                help()
                return
        if out == None:
            out = 'filter_' + str(rcut).split('.')[1]
        filterq()
    elif args[1] == '-h' or args[1] == '--h':
        help()
        return
    else:
        print "nonsensical arguments"
        help()
        return

def help():
    print "filter fastq files by their phred q-scores\n\
    \033[1mGeneral Usage\033[0m\tThis program has two modules - Qscore Evaluation, and Qscore Filtering.\n\
    \t\tYou must run Evaluation before attempting Filtering, but may filter multiple\n\
    \t\ttimes after that. \033[1mEvaluation\033[0m calculates each read's phred Qscore total sum,\n\
    \t\tand rewrites bases with individual uncertainty greater than bcut(defaults to .99)\n\
    \t\tas 'N'. If base was cut, its contribution to expected probability is bcutcomp\n\
    \t\t(defaults to .05). Evaluation outputs a log for Filter's use, and a visual\n\
    \t\tsummary of 'N' positions to a folder(defaults to {input}_{bcut}_{bcutcomp}).\n\
    \t\t\033[1mFiltering\033[0m separates out all reads with a Qscore lower than rcut, to out\n\
    \t\t(defaults to filter_{rcut}), maintaining order, with discarded sequences in\n\
    \t\tout.cut,and creates a read Qscore ofhistogram.\n\
    \033[1mEvaluation Usage\033[0m\tqScoreQC.py --e input.fastq [-f folder][-bq bcut][-bc baseCutCompensation][-nbg][nrg]\n\
    \t-f\tspecify folder evalq outputs logs to; will create folder if nonexistant\n\
    \t-bq\tspecify Qscore threshold above which individual base is replaced by 'N'\n\
    \t-bc\tspecify amount added to total read Qscore from 'N'\n\
    \t-nbg\trequest not to create graph of 'N' distribution along read positions\n\
    \t-nrg\trequest not to create histogram of number of reads with certain Qscores\n\
    \033[1mFilter Usage\033[0m\tqScoreQC.py --f folder [-o out][-rc rcut][-ng]\n\
    \t-o\tspecify name of filtered sequences output\n\
    \t-rc\tspecify total Qscore threshold above which reads not included in out"

def evalq():
    if os.path.exists(folder):
        print "warning: outputting into " + folder + ", which already exists"
    else:
        os.makedirs(folder)
    if rGraph:
        p1 = PdfPages(folder + "/qhist.pdf")
        qvals = []
        maxexp = 0
    if bGraph:
        p2 = PdfPages(folder + "/nDist.pdf")
        ndist = []
    with open(inf, 'r') as infile, open(folder + "/log", "w") as log, open(folder + "/n", "w") as myfile:       
        for line in infile:
            q = 0
            label = line
            seq = infile.next()
            plus = infile.next()
            qscore = infile.next()
            if len(seq.rstrip('\n')) != len(qscore.rstrip('\n')):
                return "sequence and qscore sequence different lengths: " + label
            for i in range(len(qscore)-1):
                score = asc2p(qscore[i])
                if score > bcut:
                    seq = list(seq)
                    seq[i] = 'N'
                    if bGraph: 
                        while len(ndist) <= i:
                            ndist.append(0)
                        ndist[i] += 1.0
                    seq = "".join(seq)
                    q += bcutcomp
                else:
                    q += score
            log.write(str(q) + "\n")
            if rGraph: 
                qvals.append(q)
                if q > maxexp:
                    maxexp = q
            myfile.write(line)
            myfile.write("".join(seq))
            myfile.write(plus)
            myfile.write(qscore)
    if rGraph:
        plot = plt.figure()
        plt.title("Reads per Cummulative Qscore\nbcut:" + str(bcut) + " comp:" + str(bcutcomp))
        plt.tick_params(axis = 'both', labelsize = 8)
        plt.ylabel("Expected Number of Errors per Read")
        plt.xlabel('Cummulative Number of Reads')
        plt.hist(qvals, 257, color = 'blue', histtype='step', orientation = 'horizontal', cumulative=True, alpha = .6)
        plt.ylim([-.01, maxexp])
        plt.gca().invert_yaxis()
        #plt.gca().invert_xaxis()
        p1.savefig(plot)
        p1.close()
    if bGraph and ndist:
        plot = plt.figure()
        plt.tick_params(axis = 'both', labelsize = 8)
        plt.ylabel("Number N / Total Reads at this Position")
        plt.xlabel('Position on Read')
        tmp = []
        for i in range(len(ndist)):
            ndist[i] /= len(qvals)
            tmp.append(i)
        plt.plot(tmp, ndist, color='blue', linestyle='dashed', marker='o',markerfacecolor='cyan', markersize=6, alpha = .6)
        plt.title("N Positions on Reads\nbcut:" + str(bcut) + " comp:" + str(bcutcomp))
        p2.savefig(plot)
        p2.close()

def filterq():
    if not os.path.exists(folder):
        print folder + " can't be found!"
        return
    with open(folder + "/n", "r") as infile, open (folder + "/log", "r") as log,\
        open(folder + "/" + out, "w") as good, open(folder + "/" + out + ".cut", "w") as bad:
        for line in log:
            if float(line) < rcut:
                good.write(infile.next())
                good.write(infile.next())
                good.write(infile.next())
                good.write(infile.next())
            else:
                bad.write(infile.next())
                bad.write(infile.next())
                bad.write(infile.next())
                bad.write(infile.next())

def asc2p(asc):
    return q2p(asc2q(asc))

def asc2q(asc):
    return ord(asc) - 33

def p2q(p):
    return -10 * math.log(p, 10)

def q2p(q):
    return 10**(-float(q) / 10)

if __name__ == '__main__':
    main()