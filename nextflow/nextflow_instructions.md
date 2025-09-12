## If using L1EM the first time:
1. Install L1EM and associated dependencies, as described in the manual.
2. Build L1EM reference, as described in the manual.

## To process your samples using the Nextflow pipeline:
3. If you already have a copy of the original L1EM scripts on your system, you just need to copy the content of the *L1EM/nextflow/* directory as well as copy the following new utility scripts:
    - *L1EM/utilities/L1EM_readpairs_nf.py*
    - *L1EM/utilities/report_l1_exp_counts_unstranded_nf.py*
    - *L1EM/utilities/report_l1_exp_counts_nf.py*
    - *L1EM/utilities/report_l1hs_transcription_unstranded_nf.py*
    - *L1EM/utilities/report_l1hs_transcription_nf.py*
4. Create csv file containing all the files you wish to process (see *L1EM/nextflow/filelist_example.csv*).
5. Prepare a directory where you'll place all your L1EM outputs. This includes:
    - Creating the directory
    - Copying *L1EM/nextflow/l1em_main.nf*, L1EM/nextflow/*l1em_config.config* and *L1EM/nextflow/run_l1em_nf_array.sh* files to that directory.
    - Add your email address to the *l1em_config.config* and *run_l1em_nf_array.sh* files at the appropriate places for notifications from slurm (search for the string "email@email.org").
    - Note: If you change the filename of *l1em_main.nf* or *l1em_config.config*, do not forget to update the filenames in the *run_l1em_nf_array.sh* file as well. 
6. Run the following commands:
    
    ```cd <directory containing l1em_main.nf, l1em_config.config, and run_l1em_nf_array.sh files>```
    
    ```sbatch --array=1-N%L run_l1em_nf_array.sh <full path to csv file containing the samples to process>```

    - N is the number of lines (excluding header) in the csv file.
    - L is the maximum number of parallel jobs you want to run. I recommend L=2 if you want to do other work on HPC in parallel, or L=5 (or more) if you want to dedicate all of your resources to L1EM. Do not run the array job without specifying L if you have >10 input files!
    

Notes:
- At the moment, nextflow jobs are submitted for each input file individually by the run_l1em_nf_array.sh script. This is not optional!!!
- If you want to make sure you're requesting the correct amount of resources for each step, I encourage you to process a couple of input files first and check the trace files in the result directory for details on how much time, memory and CPU was used by each process. Most bam files do not require the amount of resources I specified.
- The pipeline is currently set up to create index files (*.bam.bai) for the input bams, at the same location where the input bams are originally located. Please make sure that you have a copy of the original bai files, if needed.
