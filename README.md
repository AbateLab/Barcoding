# BarcodingQC
Tools for evaluating quality of barcoding runs, and filtering out low quality data.
Running any program with no parameters, or "--h" will display a general guide to its use.

qScoreQC.py:
We decided to filter on expected phred q-value instead of an average, as expected value carries more intuitive information with regards to eventual assembly quality. Expected phred q-value, or the expected number of incorrect bases per read, is calculated by summing each q-value associated with each base in a read from a fastq format file. However, because some assemblers might value an "unknown base" marking of "N", we provide the option of setting a "bcut", where bases with a q-score higher than "bcut" are replaced by "N", and because this base is no longer explicitely wrong, we provide the option of setting a "bcutcomp", which is summed into the total expected value, instead of the no longer relevant q-score.\n
qScoreQC.py --f {...} is dependant on logs calculated in qScoreQC.py --e {...}, but qScoreQC.py --e {...} needs only be run once, after which one can create different cuts of the original file with qScoreQC.py --f {...} as much as one pleases. qScoreQC.py --e {...} takes O(#reads * len(read)) time, while qScoreQC.py --f {...} takes O(#reads).\n
TODO: parallelize qScoreQC.py --e
