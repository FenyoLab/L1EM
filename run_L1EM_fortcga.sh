#!/bin/bash

# Usage: bash run_L1EM.sh /fullpathto/alignments.bam /fullpathto/L1EM /fullpathto/L1HS.fa

# Parameters
threads=16 #How many threads to use for samtools, bwa and L1EM
realignNM=3 #Number of mismatches allowed in bwa realignment
L1EM_NM=2 # Number of mismatches allowed when enumerated candidate alignments
NMdiff=1 #Skip candidate alignments with greater than this many more mismatches than the best alignment
bwa_i=20 #bwa i parameter prevents indels near the edges of a read
error_prob=0.01 #Probability of a read error at a given position
max_start2start_len=500 #Max allowed template/fragment length
reads_per_pickle=10000 #Number of rows in each G(R) matrix chunk. Decrease if memory usage is too high.
EM_threshold=1e-6 #Keep taking EM steps until no entry in X changes by more than this value. Increasing this parameter will shorten run time.
template_fraction=0.0001 #Fraction of reads to consider when calculated median template length.

# If you need to specify paths to required packages
bwa=$(which bwa) # version 0.7.17 tested
samtools=$(which samtools) # version 1.9 tested
python=$(which python) # use version 2.7

# Command line arguments
bamfile=$1
L1EM_directory=$2
L1HS=$3

L1EM_bed=$L1EM_directory'/annotation/L1EM.400.bed'
L1EM_fa=$L1EM_directory'/annotation/L1EM.400.fa'
L1EM_code_dir=$L1EM_directory'/L1EM/'
L1EM_utilities_dir=$L1EM_directory'/utilities/'

# Try to realign unaligned reads using bwa aln.
echo 'STEP 1: realign'
mkdir idL1reads
cd idL1reads
$samtools view -@ $threads -b -F 2 $bamfile | $samtools sort -@ $threads -n - | $samtools fastq - -1 unaligned.fq1 -2 unaligned.fq2
$bwa aln -k $realignNM -n $realignNM -t $threads -i $bwa_i $L1HS unaligned.fq1 > 1.sai
$bwa aln -k $realignNM -n $realignNM -t $threads -i $bwa_i $L1HS unaligned.fq2 > 2.sai
$bwa sampe $L1HS 1.sai 2.sai unaligned.fq1 unaligned.fq2 | $samtools view -b -f 2 -@ $threads - | $samtools sort -@ $threads - > realigned.bam 

# Extract L1HS/L1PA* aligning reads.
echo 'STEP 2: extract'
$python ${L1EM_utilities_dir}read_or_pair_overlap_bed.py $L1EM_bed $bamfile temp.bam
$samtools sort -@ $threads -n temp.bam | $samtools fastq - -1 L1.fq1 -2 L1.fq2
$samtools fastq realigned.bam -1 temp.fq1 -2 temp.fq2
cat temp.fq1 >> L1.fq1
cat temp.fq2 >> L1.fq2
rm temp*

# Split the L1 fastq files for parallel execution
cd ..
mkdir split_fqs
split_fq_size=$(wc -l idL1reads/L1.fq1 | awk '{print $1/('$threads'*4)+1}' | cut -d '.' -f 1 | awk '{print $1*4}')
split -l $split_fq_size idL1reads/L1.fq1 split_fqs/L1.fq1.
split -l $split_fq_size idL1reads/L1.fq2 split_fqs/L1.fq2.
cd split_fqs

# Generate candidate alignments
echo 'STEP 3: candidate alignments'
for name in *.fq1.*
    do reads1=$name
    reads2=$(echo $name|sed 's/fq1/fq2/g')
    ref=$L1EM_fa
    base=$(echo $name|sed 's/.fq1//g')
    $bwa aln -t $threads -N -n $L1EM_NM -k $L1EM_NM -i $bwa_i -R 10000000 $ref $reads1 > $base.R1.aln.sai
    $bwa aln -t $threads -N -n $L1EM_NM -k $L1EM_NM -i $bwa_i -R 10000000 $ref $reads2 >  $base.R2.aln.sai
done
for name in *.fq1.*
    do reads1=$name
    reads2=$(echo $name|sed 's/fq1/fq2/g')
    ref=$L1EM_fa
    base=$(echo $name|sed 's/.fq1//g')
    $bwa sampe -n 10000000 -N 10000000 $ref $base.R1.aln.sai  $base.R2.aln.sai $reads1 $reads2 | $samtools view -bS - | $samtools sort -n - > $base.aln.bam &
done
wait

# Make G_of_R matrix
echo 'STEP 4: G(R) matrix construction'
mkdir ../G_of_R
cd ../G_of_R
medianinsert=$($python ${L1EM_utilities_dir}median_template.py $bamfile $template_fraction)
for bam in ../split_fqs/*.bam
	do $python ${L1EM_code_dir}G_of_R_unstranded.py -b $bam -i $medianinsert -p $(echo $bam| cut -d '/' -f 3) -e $error_prob -m $max_start2start_len -r $reads_per_pickle -n $NMdiff &
done
wait

# RUN EM
echo 'STEP 5: Expectation maximization'
mkdir ../L1EM/
cd ../L1EM/
ls ../G_of_R/*pk2 > G_of_R_list.txt
cp $(ls ../G_of_R/*TE_list.txt | head -1) TE_list.txt
python ${L1EM_code_dir}L1EM.py -g G_of_R_list.txt -l TE_list.txt -t $threads -s $EM_threshold

#Write results as text file
echo 'STEP 6: Writing results'

$python ${L1EM_utilities_dir}report_l1_exp_counts_unstranded.py > ../full_counts.txt
$python ${L1EM_utilities_dir}report_l1hs_transcription.py > ../l1hs_transcript_counts.txt

#Clean up
echo 'STEP 7: Clean up'
cp *final.pkl ../
cd ..

#rm idL1reads/*
#rmdir idL1reads
#rm split_fqs/*
#rmdir split_fqs
#rm G_of_R/*
#rmdir G_of_R
#rm L1EM/*
#rmdir L1EM

