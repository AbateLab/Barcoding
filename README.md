# Barcoding
Tools for evaluating quality, and filtering barcoding data. Still in Development

Run any script without arguments, or -h or --help for basic usage description.
###qScoreQC.py:
qScoreQC provides several sequencing read quality control options:

1. Replace bases below a certain q-score with 'N', and replace the associated q-score with specified 'N' q-score value
2. Cut the end of the read if the read is actually short, and Illumina just put N's to fill in spaces (still in the rough)
3. Filter reads using "Expected Number of Incorrect Bases per Read", or the sum of the probabilities of each base in a read being incorrect. Reads with Exp < specified amount(0 -> len(read)) will be placed in the good pile.
4. Filter reads using "Minimum Allowable Q-Score", where only reads with all individual-base q-scores > specified q-score will be placed in the good pile.
5. Filter reads using "Minimum Percent of Read > specified Q-Score", where only reads with percent of reads > specified Q-Score > specified percent will be placed in the good pile.

Notice several dependencies between these options which influence workflow - First, all the filter options depend on what bases are replaced with N's, which modifies q-scores. This is only reasonable, because the information carried by "N" is less precise, but more accurate, and Q-score reflects accuracy. Second, options 2 through 5 may be repeated multiple times with cutoff changes, but their heuristic scores would not change between these repetitions.

In light of this, we decide to accomplish qScoreQC.py's objective using two modules:

1. (Eval) Replace all the N's needed, Calculate all the Heuristics. This depends on a folder to save logs into, an input file, an N cutoff, an N replacement Q-Score, and the specified Q-Score for option 5. (Only need to do once)
2. (Filter) Filter/trim reads in whatever way you wish. This depends on a folder with logs, an output tag, and whatever filter/trim cutoff you choose.
Eval is O(read length * reads), and filter is O(reads)

To assist with choosing the right parameters for filtering, replacing, and trimming, qScoreQC provides a few graphs, also saved into the folder specified in the eval run. After the eval run, five graphs can be provided with the '-g' tag - a cummulative histogram of number of N's a read ends with (will be changed as 2. changes), a base composition per read position graph, and cummulative histograms of scores for each heuristic. On each filter run, a histogram of the scores remaining, using the filter run's selected heuristic, will be provided. Also, general reports can be found in folder/report.

###dfsCluster.py:
Intended to cluster pcr output sequences into groups of sequences all originating from a unique pcr input sequence

As the name suggests, dfsCluster employs a dfs over a hamming space, but can only visit nodes specified by pcr output sequences. This should work because

1. PCR in the overwhelming majority of the time will not make more than 1 error per sequence
2. The hamming space is large enough, and the error rate is low enough, that random walks originating from random points in the space will not ever come closer than 2 hamming distance of another path.

Currently creating simulation barcodes to confirm this is valid

dfsCluster reads a fasta/q file of reads (nonunique), groups reads into barcode groups (unique), removes barcode groups with fewer reads than a specified value (default 1), and groups the barcodes into larger clusters assigned to a "center", or the sequence composed of the majority base at each position in the barcodes in the cluster. dfsCluster then creates a jackpottogram (pie chart with each slice being the number of reads in a single cluster), a histogram and pie chart binned by cluster size, a histogram binned by minimum hamming distance between a center and all other centers, and a .cid file which prints all the id's of reads associated with a cluster below the cluster's center. Ex:

>\>"center1"  
[id1, id2, id3, ...]  
\>"center2"  
[id10, id11, ...]

It is recommended to run with -v (verbose), or at least save the terminal output to a file, because information such as #unique barcodes, #barcodes cut, %reads cut, #clusters formed will be printed.

TODO:Everything

###jackpotQC.py:
jackpotQC provides visualizations for "jackpotting" of barcode groups, or the quality of PCR amplifying different sequences in hugely uneven quantities. Input a fasta/q file of sequences, and jackpottogram will create a histogram, and a pie chart based on the number of copies of each unique sequence.

TODO:Make cummulative histogram binned by barcode size

###readLenHist.py:
Just a simple histogram of read length maker. Nothing to it. Tell it what fastq file to graph and the name of the output pdf.

TODO:Open to suggestions
