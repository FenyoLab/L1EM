# Parameters
export threads=16 #How many threads to use for samtools, bwa and L1EM
export realignNM=3 #Number of mismatches allowed in bwa realignment
export L1EM_NM=3 # Number of mismatches allowed when enumerated candidate alignments
export NMdiff=2 #Skip candidate alignments with greater than this many more mismatches than the best alignment
export bwa_i=20 #bwa i parameter prevents indels near the edges of a read
export error_prob=0.01 #Probability of a read error at a given position
export max_start2start_len=500 #Max allowed template/fragment length
export reads_per_pickle=10000 #Number of rows in each G(R) matrix chunk. Decrease if memory usage is too high.
export EM_threshold=1e-7 #Keep taking EM steps until no entry in X changes by more than this value. Increasing this parameter will shorten run time.
export template_fraction=1 #Fraction of reads to consider when calculated median template length.
