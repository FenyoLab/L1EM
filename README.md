## Installation
### conda way
You will need
1. git (https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
2. anaconda (https://docs.anaconda.com/anaconda/install/)

Download from github
```
git clone https://github.com/FenyoLab/L1EM
```
Create conda environment
```
cd L1EM
conda env create -f L1EM.yml
```

Before running L1EM, activate the environment:
```
source activate L1EM
```

When finished, deactivate the environment:
```
source deactivate L1EM
```

### old way
Alternatively you can install the following dependencies yourself:
* python version 2.7+ (version 2.7 tested)
* bwa (version 0.7.17 tested)
* samtools (version 1.9 tested)
* numpy (version 1.14.3 tested)
* scipy (version 1.1.0 tested)
* pysam (version 0.15.0 tested)
* bedtools (version 2.27.1 tested)

No compiling of L1EM is necessary. Python scripts will be called from inside the L1EM
directory.

If necessary, you can specify the path for bwa and samtools in the run\_L1EM.sh script.
You must use samtools >=1.0. Early version of pysam will not work. I highly recommend
that you use bwa 0.7.17. Earlier versions may differ in how they write the XA tag. This
will lead to inaccurate results without throwing an error.

## Quick guide
### First time: build L1EM reference
You will need the hg38 reference genome in fasta format, with bwa index.
Downloaded from UCSC genome browser:
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
zcat hg38.fa.gz > hg38.fa
bwa index hg38.fa
```
Note: this will take some time.

Then you can build the L1EM reference using the provided shell script:
```
bash generate_L1EM_fasta_and_index.sh /fullpathto/hg38.fa
```
This should be done inside the L1EM directory

### Executing the L1-EM pipeline
You will need a bam file with strand specific paired end read alignments to hg38. You can
use any aligner, but make sure that all reads from the original fastq files are present
trimming should be okay, but is tested. Filtering reads will potentially break the pipeline.

First move to an empty directory and then execute the shell script:
```
bash -e /fullpathto/run_L1EM.sh /fullpathto/alignments.bam /fullpathto/L1EM /fullpathto/hg38.fa
```
L1EM will write files with specific names, so do NOT run two instances of L1EM in the same
directory.

At the end of the run\_L1EM.sh script are a commented set of commands to delete all the
intermediate files. If you wish to automatically delete intermediate files, you can delete
these comments.

### Output
At completion, three tab delimited tables will be written.
1. full\_counts.txt: raw count estimates for each L1HS/L1PA\* element with any aligned read pairs
2. l1hs\_transcript\_counts.txt: expression estimates for L1HS elements, reported as raw counts
3. filter\_L1HS\_FPM.txt: L1HS whose expression is supported by at least 100 read pairs, reported as FPM (read pairs per million properly aligned)

The rows of all files are L1 loci.

For full\_counts.txt each of the five transcript types:
only, runon, passive (sense), passive (antisense), antisense
are reported.

For l1hs\_transcript\_counts.txt and filter\_L1HS\_FPM.txt only proper transcription from L1HS elements start at the
5' UTR is reported.

The results are also written as pickle files to facilitate further analysis in python. To
generate a python dictionary with keys being the transcript names and values being the
relative expression:
```
X_est = dict(zip(pickle.load(open('names_final.pkl')),pickle.load(open('X_final.pkl'))))
```

## Additional details
* Our Bioinformatics paper introducing L1EM: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz724/5581349
* More details can be found in manual.md

## Mouse Version
Scripts and annotation to measure the expression of LINE-1 loci in mm39 has been added. The mouse version uses all the same methodology as the human version, but has not been as rigorously tested.
1. Download and index the mm39 reference genome (UCSC genome browser version)
```
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
zcat mm39.fa.gz > mm39.fa
bwa index mm39.fa
```
2. Build the mm39 L1EM reference.
```
bash generate_mm39_L1EM_fasta_and_index.sh /fullpathto/mm39.fa
```
3. Run L1EM.
```
bash /fullpathto/run_L1EM_mm39.sh /fullpathto/alignments.bam /fullpathto/L1EM /fullpathto/mm39.fa
```
All L1Md loci are quantified in full\_counts.txt. Normalized expression of 5' UTR intact young (L1Md\_Tf I/II/II, L1Md\_Gf I/II, L1Md\_A I/II/III) LINE-1 loci supported by at least 100 reads can be found in filter\_active\_L1Md\_FPM.txt.




