## Pipeline Parameters

The key parameters for L1EM are listed at the beginning of the run\_L1EM.sh file. Default parameters should work well in most cases, but advanced users may wish to tinker.
1. threads. Dictates the number of threads that L1EM will spawn. More threads will improve parallel performance, but memory usage scales linearly with number of threads.
2. realignNM. The number of mismatches to allow when trying to realign reads that do not align as proper pairs in the bam file provided. Default is 3, but you might want to increase for longer reads.
3. L1EM_NM. As above, but for the generation of candidate alignments to the L1EM reference. Including more candidate alignments will slow the computation, but too few candidate alignments could yield less accurate results.
4. NMdiff. Only consider alignments with at most this many more mismatches than the primary alignment. Because read likelihood diminishes exponentially with additional mismatches, increasing this parameter is unlikely to affect results but will slow the EM steps.
5. bwa\_i. By default bwa will create a large number of alignments with indels near the edge of the read. This parameter will prevent this behavior. You may wish to decrease this parameter for shorter reads.
6. error\_prob. Probability of an error. Error probability is chosen to be constant because computing the read likelihood from base quality scores is slow.
7. max\_start2start\_len=500. Maximum allowed fragment/template length. Increase if you are using data with very large fragments.
8. reads\_per\_pickle. The G(R) matrix is split into a number of pickle files, so the entire matrix doesn't need to sit in memory. Decreasing this parameter will free up memory at the G(R) construction and EM steps.
9. EM\_threshold. Run EM steps until no entry in X changes by more than this value. The paremeter is chosen to be small by default to ensure convergence. Increasing the parameter modestly will improve run time.
10. template\_fraction. When computing median template length, subsample read to this fraction. You only need about 10,000 proper pairs to get a good estimate.

## Generating new annotations
If you wish to run L1-EM for another retrotransposon or for another model organism, you will need to generate a new annotation.
1. Create a bedfile with the following naming scheme:
family.category.region.strand
Where family is the name of the repeat family,
category is 1 is the element has a promoter and 0 otherwise
region is the genome region (chrom:start-stop) of the element
strand is +/- depending which strand the element falls on
The bedfile must have the six required fields: chrom, start, stop, name, score, strand
The start and stop coordinates should include 400 positions of flanking sequence on either end.
Exons overlapping the annotation can also be included.
2. Create a fasta file from your bed file and index it with bwa:
```
bedtools getfasta -s -name -fi refernece.fa -bed annotation.bed > annotation.fa
bwa index annotation.fa
```
3. Update lines 27 and 28 to point toward your new annotation.

## Pipeline steps
### STEP 1: realign
In this step reads that are not properly paired are extracted and realigned with bwa. Many aligners do not bother with highly redundant reads, so this step is included to ensure that LINE-1 aligning reads are identified.

### STEP 2: extract
In this step, L1HS/L1PA reads are extracted. Any read pair for which either end overlaps an entry in the L1EM.400.bed annotation is considered.

### STEP 3: candidate alignments
The extracted reads are aligned to L1EM.400.fa, all secondary alignments with up to L1EM_NM mismatches are found. The candidate alignments fastqs are split for parallelization. It is vitally important that all candidate alignments are identified. Missing some of these alignments will drastically hurt accuracy. For this reads bwa aln is used. Do not use bwa mem or STAR as these aligner do not provide a complete enumeration of secondary alignments for highly repetitive elements (like L1HS).

### STEP 4: G(R) matrix construction

The bam files of candidate alignments are read by the script G\_of\_R.py. The likelihood of each candidate alignment is calculated and added to the G(R) matrix.

The following options are additional parameters that can be accessed at this step:
1. -f/--flanking specifies the amount of flanking sequence in the annotation. If you created you own annotation with more or less that 400 bases of flanking sequence specify that here.
2. --as\_start. If you wish to change to TSS for antisense transcription do that here.
3. -w/--wiggle. Some proper LINE-1 transcripts start slightly before the annotation start of the 5'UTR. This parameter extends the annotated element this many bases in either direct (default is 20).
4. --min\_len. Puts a floor on transcript effective length to prevent cases where transcription of very short elements are over predicted. Default is 500.
5. --min\_exon\_len. Corresponding minimun effective length for exon annotations. Default is 100.

### STEP 5: Expectation maximization
In this step, the expectation maximization algorithm is used to compute a maximum likelihood estimate of relative expression, using the G(R) matrix output in the previous as input.
The following options are additional parameters that can be accessed at this step:
1. -r/--report\_every. Write the estimate every n steps.
2. -m/--max\_nEMsteps. By default EM stops if converge has not been achieved after 10000 steps. Change that value here.

### STEP 6: Writing results
At completion, two tab delimited table will be written.
1. full\_counts.txt
2. l1hs\_transcript_counts.txt

### STEP 7: Clean up
All the intermediate files are delete at this step. Comment out these lines if you want to keep them.

The rows of both files are L1 loci.

For full\_counts.txt each of the five transcript types:
only, runon, passive (sense), passive (antisense), antisense
are reported.

For l1hs\_transcript_counts.txt only proper transcription from L1HS elements start at the
5' UTR is reported.

The results are also written as pickle files to facilitate further analysis in python. To
generate a python dictionary with keys being the transcript names and values being the
relative expression:
```
X_est = dict(zip(pickle.load(open('names_final.pkl')),pickle.load(open('X_final.pkl'))))
```




