#!/bin/bash
#SBATCH --job-name=nf_job
#SBATCH --output=logs/l1em_nf_%j.out
#SBATCH --mail-user=email@email.org
#SBATCH --mail-type=ALL
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu_medium,cpu_short,cpu_long

# Usage: sbatch --array=1-N%L run_l1em_nf_array.sh <full path to csv file>
# where N is the number of lines (excluding header) in the csv file. See example: filelist_example.csv
# L is the maximum number of concurrent jobs you want to run
# (I recommend L=2 if you want to do other work on HPC in parallel, or L=5 if you want to dedicate all of your resources to L1EM)
# Do not forget to add your email address to the --mail-user field above to get notifications from slurm!

# TO DO: Right now we need to enter a specific directory containing the bash and nextflow files to run this script.
# In the future, think about ways to make this more flexible (e.g. passing the directory as an argument to the script)
# Also, if the filenames of l1em_main.nf or l1em_config.config are changed, the script needs to be updated accordingly

# TO DO: after the array finished, think about ways to automatically check if all jobs were successful -- if not,
# provide their information in a table for review (why did it fail; can we just resubmit it or do we have to change something?)

module load nextflow/24.10.3 
module load python/cpu/2.7.15-ES
module load bwa/0.7.17
module load samtools/1.9
module load bedtools/2.27.1

PATH="/usr/local/bin:$PATH"

SAMPLESHEET="$1"

# Get the line for this array task (skip header)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLESHEET")

# Read columns into variables
IFS=, read -r sample_id strandedness input_bam hg38_fa_dir ref l1em_scripts results_dir nr_threads realignNM L1EM_NM NMdiff bwa_i error_prob max_start2start_len reads_per_pickle EM_threshold template_fraction <<< "$LINE"

mkdir -p $results_dir
chmod 770 $results_dir
cp l1em_config.config "$results_dir"/
cp l1em_main.nf "$results_dir"/
# TO DO: check if this there is a better way to process files individually without copying all of these files each time 

echo "Running Nextflow pipeline for sample: $sample_id"

cd $results_dir

nextflow -C "$results_dir"/l1em_config.config run "$results_dir"/l1em_main.nf -w "$results_dir"/work_tmp -with-trace "$results_dir"/trace_$(date +%Y%m%d-%H%M%S).txt \
  --sample_id "$sample_id" \
  --strandedness "$strandedness" \
  --input_bam "$input_bam" \
  --hg38_fa_dir "$hg38_fa_dir" \
  --ref "$ref" \
  --l1em_scripts "$l1em_scripts" \
  --results_dir "$results_dir" \
  --nr_threads "$nr_threads" \
  --realignNM "$realignNM" \
  --L1EM_NM "$L1EM_NM" \
  --NMdiff "$NMdiff" \
  --bwa_i "$bwa_i" \
  --error_prob "$error_prob" \
  --max_start2start_len "$max_start2start_len" \
  --reads_per_pickle "$reads_per_pickle" \
  --EM_threshold "$EM_threshold" \
  --template_fraction "$template_fraction"


# If the nextflow run was successful, delete the work_tmp directory
if [ $? -eq 0 ]; then
  echo "Nextflow run completed successfully."
  rm -r "$results_dir"/work_tmp

else
  echo "Nextflow run failed."
fi


