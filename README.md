# Barcoding
Tools for evaluating quality, and filtering barcoding data.

###qScoreQC.py:
ASSUMES FASTQ and FASTA "reads" and "quality" are single lines!

We judge each read by Expected Number of Incorrect Bases, or the sum of each base's probability of being incorrect taken fron the fastq file's phred q-score.

We provide an option to cut a base after a certain uncertaintly level, and replace it with an N, and the user can set a "probability incorrect" value to add to the Expected Number of Incorrect Bases score.

We provide an option to cut all sequential N's off the end of a read if there are more than a user specified value. Please note that most existing bioinformatics tools don't play nice with files with non-homogenous read lengths.

We provide the capability of graphing certain characteristics of read quality - a cummulative histogram of read scores, a graph of the number of reads with score below a user set threshold, a graph of the distribution of N's over the positions on the read, and a cummulative histogram of reads ending in x number of sequential N's

One must run the eval option before other options are possible. One only needs to run the eval option once. Please read qScoreQC's help text for more specific usage help.

TODO: expand "cut the end of the read if the read is actually short, and Illumina just put N's to fill in spaces" functionality.

###dfsCluster.py:
Intended to cluster pcr output sequences into groups of sequences all originating from a unique pcr input sequence

Still in the works as to UI, although it will print in terminal the clusters it finds. As the name suggests, dfsCluster employs a dfs over a hamming space, but can only visit nodes specified by pcr output sequences. This should work because

1. PCR in the overwhelming majority of the time will not make more than 1 error per sequence
2. The hamming space is large enough, and the error rate is low enough, that random walks originating from random points in the space will not ever come closer than 2 hamming distance of another path.

TODO:Everything

###readLenHist.py:
Just a simple histogram of read length maker. Nothing to it. Tell it what fastq file to graph and the name of the output pdf.

TODO:Open to suggestions
