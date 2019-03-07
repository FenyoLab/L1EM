## Installation
### conda way
You will need
1. git (https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
2. anaconda (https://docs.anaconda.com/anaconda/install/)

Download from github
```
git clone https://github.com/wmckerrow/L1EM
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
* python version 2.7
* bwa (version 0.7.17 tested)
* samtools (version 1.9 tested)
* numpy (version 1.14.3 tested)
* scipy (version 1.1.0 tested)
* pysam (version 0.15.0 tested)
* bedtools (version 2.27.1 tested)

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

### Executing the L1-EM pipeline
You will need a bam file with strand specific paired end read alignments to hg38. You can
use any aligner, but make sure that all reads from the original fastq files are present
trimming is okay, but filtering reads will potentially break the pipeline.

First move to an empty directory and then execute the shell script:
```
bash /fullpathto/run_L1EM.sh /fullpathto/alignments.bam /fullpathto/L1EM /fullpathto/hg38.fa
```
L1EM will write files with specific names, so do NOT run two instances of L1EM in the same
directory.

At the end of the run\_L1EM.sh script are a commented set of commands to delete all the
intermediate files. If you wish to automatically delete intermediate files, you can delete
these comments.

### Output
At completion, two tab delimited table will be written.
1. full\_counts.txt
2. l1hs\_transcript_counts.txt

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

## Additional details
* URL for our methods paper will appear here when it is published
* More details can be found in manual.md



