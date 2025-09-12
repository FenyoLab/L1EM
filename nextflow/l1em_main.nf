#!/usr/bin/env nextflow


Channel
  .of(
    tuple(
      params.sample_id, // sample name
      params.strandedness, // "stranded" or "unstranded"
      file(params.input_bam), // input BAM file
      file(params.hg38_fa_dir), // directory containing hg38.fa and its index files
      file(params.ref), // reference genome, as prepared based on L1EM instructions
      file(params.l1em_scripts), // directory containing L1EM scripts
      params.results_dir, // directory to store final results
      params.nr_threads.toInteger(), // number of threads to use
      params.realignNM.toInteger(), // realignNM parameter, as specified in the L1EM documentation
      params.L1EM_NM.toInteger(), // L1EM_NM parameter, as specified in the L1EM documentation
      params.NMdiff.toInteger(), // NM_diff parameter, as specified in the L1EM documentation
      params.bwa_i.toInteger(), // bwa_i parameter, as specified in the L1EM documentation
      params.error_prob, // error_prob parameter, as specified in the L1EM documentation
      params.max_start2start_len.toInteger(), // max_start2start_len parameter, as specified in the L1EM documentation
      params.reads_per_pickle.toInteger(), // reads_per_pickle parameter, as specified in the L1EM documentation
      params.EM_threshold, // EM_threshold parameter, as specified in the L1EM documentation
      params.template_fraction // template_fraction parameter, as specified in the L1EM documentation
    )
  )
  .set { samples_ch }
// TO DO: add sanity checks for input files and parameters
// TO DO: some of these parameters are redundant (hg38_fa_dir and ref). Could be simplified in the future



workflow {
  // TO DO: better way to specify the parameters. Right now, everything is explicitly passed through the channels, which clutters the code a bit
  // TO DO: edit processes so that multiple input files can be handled (right now this is not correctly set up in combine_gofr_outputs)

  // TO DO: add option to skip setup if the input BAM is already indexed OR make sure that the index file ends up in the temporary workdir of nextflow instead of the location of the original bam file
  setup_done_ch = samples_ch
    .map { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
           nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
           max_start2start_len, reads_per_pickle, EM_threshold, template_fraction ->
           
      tuple(sample_id, input_bam)
    } | setup



    realign_input_ch = samples_ch
        .join(setup_done_ch)
        .map { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
             nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
             max_start2start_len, reads_per_pickle, EM_threshold, template_fraction,setup_done, bai_file ->

            tuple(sample_id, setup_done, input_bam, bai_file, nr_threads, realignNM, bwa_i, hg38_fa_dir)
        }
    
    
    realign_out_ch = realign_reads(realign_input_ch)

    extract_l1_input_ch = samples_ch
        .join(realign_out_ch)
        .join(setup_done_ch)
        .map { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
             nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
             max_start2start_len, reads_per_pickle, EM_threshold, template_fraction,realigned_bam, 
             realigned_bam_bai, realignment_done, setup_done, bai_file ->

            tuple(sample_id, input_bam, l1em_scripts, template_fraction,realigned_bam, realigned_bam_bai, realignment_done, bai_file)
        }

    l1_out_ch = extract_l1_reads(extract_l1_input_ch)


    // TO DO: check if there's a better way to set up this arrayjob
    align_input_ch = samples_ch
      .join(l1_out_ch)
      .flatMap { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
            nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
            max_start2start_len, reads_per_pickle, EM_threshold, template_fraction,
            l1_fq1_list, l1_fq2_list, extract_done, baminfo ->

        // Pair up fq1 and fq2 files by index and emit a tuple for each pair
        def fq1_files = l1_fq1_list instanceof List ? l1_fq1_list : [l1_fq1_list]
        def fq2_files = l1_fq2_list instanceof List ? l1_fq2_list : [l1_fq2_list]
        fq1_files.indices.collect { idx ->
          tuple(sample_id, strandedness, l1em_scripts, fq1_files[idx], fq2_files[idx], ref, realignNM, bwa_i, nr_threads, baminfo, error_prob,
                max_start2start_len, reads_per_pickle, NMdiff)
        }
      }


    align_out = generate_alignments_and_gofr(align_input_ch)
        .map { sample_id, bam, pk2, te, done ->
            tuple(sample_id, bam, pk2, te)
        }
        .groupTuple()
        // Result: [ sample_id, [bam1, bam2, ...], [pk21, pk22, ...], [te1, te2, ...] ]
        .map { sample_id, bam_list, pk2_nested, te_list ->
            // Flatten pk2_nested so it's just a flat list of paths
            def pk2_list = pk2_nested.flatten()
            tuple(sample_id, bam_list, pk2_list, te_list)
    }


    combine_gofr_input_ch = samples_ch
        .join(align_out)
        .map { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
            nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
            max_start2start_len, reads_per_pickle, EM_threshold, template_fraction, bam_list, pk2_list, te_list  ->

            tuple(sample_id, bam_list, pk2_list, te_list )
      }
    
    combined_gofr_out = combine_gofr_outputs(combine_gofr_input_ch)


    em_input_ch = samples_ch
        .join(align_out)
        .join(combined_gofr_out)
        .map { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
             nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
            max_start2start_len, reads_per_pickle, EM_threshold, template_fraction,
            bam_list, pk2_list, te_list, 
            gofr_done, g_r_list_txt, te_list_txt ->

            tuple(sample_id, l1em_scripts, gofr_done, g_r_list_txt, te_list_txt,
                  pk2_list, te_list, EM_threshold, nr_threads)
        }

    em_out = run_em(em_input_ch)

    finalize_input_ch = samples_ch
        .join(em_out)
        .join(align_out)
        .join(combined_gofr_out)
        .join(l1_out_ch)
        .map { sample_id, strandedness, input_bam, hg38_fa_dir, ref, l1em_scripts, results_dir,
             nr_threads, realignNM, L1EM_NM, NMdiff, bwa_i, error_prob,
            max_start2start_len, reads_per_pickle, EM_threshold, template_fraction,
            em_done, final_names_pkl, final_x_pkl, bam_list, pk2_list, te_list, 
            gofr_done, g_r_list_txt, te_list_txt,
            l1_fq1_list, l1_fq2_list, extract_done, baminfo  ->

            tuple(sample_id, strandedness, l1em_scripts, results_dir, em_done,
                  baminfo, g_r_list_txt, pk2_list, final_names_pkl, final_x_pkl)
        }
    
    finalize_results(finalize_input_ch)

}


// TO DO: check if .done files are really necessary, or if we can just use the existence of the output files to indicate completion of a step
// TO DO: write tests for each process & for each python script
// TO DO: add error handling in each process (e.g. check if expected output files were created, if intermediate steps were successful, etc.)
// TO DO: check if any bash commands could be simplified and (later on) made faster


process setup {
  tag { sample }

  input:
    tuple val(sample_id), path(input_bam)

  output:
    tuple val(sample_id), path("${sample_id}_setup.done"), path("${input_bam}.bai"), emit: setup_done

  script:
  """
  samtools index ${input_bam}
  touch "${sample_id}_setup.done"
  """
}


process realign_reads {
  tag { sample }

  input:
    tuple val(sample_id), path(setup_done), path(input_bam), path(bai_file), val(nr_threads), val(realignNM), val(bwa_i), path(hg38_fa_dir)

  output:
    tuple val(sample_id), path("${sample_id}_realigned.bam"), path("${sample_id}_realigned.bam.bai"), path("${sample_id}_realignment.done"), emit: realignment_done

  script:
  """
  samtools view -@ ${nr_threads} -b -F 2 ${input_bam} | \
    samtools sort -@ ${nr_threads} -n - | \
    samtools fastq - -1 ${sample_id}_unaligned.fq1 -2 ${sample_id}_unaligned.fq2

  bwa aln -k ${realignNM} -n ${realignNM} -t ${nr_threads} -i ${bwa_i} ${hg38_fa_dir}/hg38.fa ${sample_id}_unaligned.fq1 > ${sample_id}_1.sai
  bwa aln -k ${realignNM} -n ${realignNM} -t ${nr_threads} -i ${bwa_i} ${hg38_fa_dir}/hg38.fa ${sample_id}_unaligned.fq2 > ${sample_id}_2.sai

  bwa sampe ${hg38_fa_dir}/hg38.fa ${sample_id}_1.sai ${sample_id}_2.sai ${sample_id}_unaligned.fq1 ${sample_id}_unaligned.fq2 | \
    samtools view -b -@ ${nr_threads} - | \
    samtools sort -@ ${nr_threads} -o ${sample_id}_realigned.bam

  samtools index ${sample_id}_realigned.bam
  touch "${sample_id}_realignment.done"
  """
}



process extract_l1_reads {
  tag { sample }
  input:
    tuple val(sample_id),path(input_bam),path(l1em_scripts),val(template_fraction),path(realigned_bam),path(realigned_index),path(realignment_done),path(bai_file)

  output:
    tuple val(sample_id), path("${sample_id}_L1.fq1_*"), path("${sample_id}_L1.fq2_*"), path("${sample_id}_extract.done"), path("${sample_id}_baminfo.txt"),  emit: l1_fq1

  script:
  """
  python ${l1em_scripts}/utilities/read_or_pair_overlap_bed.py ${l1em_scripts}/annotation/L1EM.400.bed ${sample_id}_realigned.bam ${sample_id}_temp.bam
  samtools sort -n ${sample_id}_temp.bam | samtools fastq - -1 ${sample_id}_temp1_L1.fq1 -2 ${sample_id}_temp1_L1.fq2
  python ${l1em_scripts}/utilities/read_or_pair_overlap_bed.py ${l1em_scripts}/annotation/L1EM.400.bed ${input_bam} ${sample_id}_temp.bam
  samtools sort -n ${sample_id}_temp.bam | samtools fastq - -1 ${sample_id}_temp2_L1.fq1 -2 ${sample_id}_temp2_L1.fq2
  cat ${sample_id}_temp1_L1.fq1 ${sample_id}_temp2_L1.fq1 > ${sample_id}_L1.fq1
  cat ${sample_id}_temp1_L1.fq2 ${sample_id}_temp2_L1.fq2 > ${sample_id}_L1.fq2
  rm ${sample_id}_temp1_L1.fq1 ${sample_id}_temp2_L1.fq1
  rm ${sample_id}_temp1_L1.fq2 ${sample_id}_temp2_L1.fq2
  split_fq_size=\$(wc -l ${sample_id}_L1.fq1 | awk '{print \$1/(16*4)+1}' | cut -d '.' -f 1 | awk '{print \$1*4}') # always split up to 16 files
  split -l \$split_fq_size ${sample_id}_L1.fq1 ${sample_id}_L1.fq1_
  split -l \$split_fq_size ${sample_id}_L1.fq2 ${sample_id}_L1.fq2_
  python ${l1em_scripts}/CGC/median_template_and_pairs.py ${input_bam} ${template_fraction} > ${sample_id}_baminfo.txt #original l1em.sh script uses the constant 0.001!!!
  touch "${sample_id}_extract.done"
  """
}


process generate_alignments_and_gofr {
  tag { sample }
  input:
    tuple val(sample_id),val(strandedness),path(l1em_scripts),path(l1_fq1),path(l1_fq2),val(ref),val(realignNM),val(bwa_i),val(nr_threads),path(baminfo),val(error_prob),val(max_start2start_len),val(reads_per_pickle),val(NMdiff)

  output:
    tuple val(sample_id),path("${sample_id}_L1_*.aln.bam"), path("${sample_id}_L1_*.aln.bam.*.pk2"), path("${sample_id}_L1_*.aln.bam_TE_list.txt"), path("${sample_id}_array_alignment.done"),  emit: l1_alignments // 

  script:
  """
  reads1=${l1_fq1}
  reads2=${l1_fq2}
  base=\$(echo ${l1_fq1}|sed 's/.fq1//g')
  fq1_base=\$(basename ${l1_fq1} | sed 's/\\.fq1.*//')
  fq2_base=\$(basename ${l1_fq2} | sed 's/\\.fq2.*//')
  if [[ "\$fq1_base" != "\$fq2_base" ]]; then
    echo "ERROR: FQ1 and FQ2 file base names do not match: \$fq1_base vs \$fq2_base" >&2
    exit 1
  fi

  bwa aln -t ${nr_threads} -N -n ${realignNM} -k ${realignNM} -i ${bwa_i} -R 10000000 ${ref} \$reads1 > \$base.R1.aln.sai
  bwa aln -t ${nr_threads} -N -n ${realignNM} -k ${realignNM} -i ${bwa_i} -R 10000000 ${ref} \$reads2 > \$base.R2.aln.sai
  bwa sampe -n 10000000 -N 10000000 ${ref} \$base.R1.aln.sai \$base.R2.aln.sai \$reads1 \$reads2 > temp.\$base.aln.sam
  samtools view -@ ${nr_threads} -bS temp.\$base.aln.sam > temp.\$base.aln.bam
  samtools sort -@ ${nr_threads} -n temp.\$base.aln.bam -o \$base.aln.bam
  rm temp.\$base.aln.sam temp.\$base.aln.bam \$base.R1.aln.sai \$base.R2.aln.sai

  medianinsert=\$(head -1 ${sample_id}_baminfo.txt)
  if [[ "$strandedness" == "stranded" ]]; then
    python ${l1em_scripts}/L1EM/G_of_R.py -b \$base.aln.bam -i \$medianinsert -p \$base.aln.bam -e ${error_prob} -m ${max_start2start_len} -r ${reads_per_pickle} -n ${NMdiff}
  else
    python ${l1em_scripts}/L1EM/G_of_R_unstranded.py -b \$base.aln.bam -i \$medianinsert -p \$base.aln.bam -e ${error_prob} -m ${max_start2start_len} -r ${reads_per_pickle} -n ${NMdiff}
  fi

  #echo "Current directory contents:"
  #ls -lh .

  touch "${sample_id}_array_alignment.done"
  """
}




process combine_gofr_outputs {
  //maxForks 1  // limits to 1 concurrent execution
  tag { sample }

  input:
    tuple val(sample_id), path(bam_list), path(pk2_list), path(te_list)
  
  output:
    tuple val(sample_id),path("${sample_id}_gofr.done"), path("${sample_id}_G_of_R_list.txt"), path("${sample_id}_TE_list.txt"), emit: combined_gofr

  script:
  """
  ls ${sample_id}_L1_*.aln.bam.*.pk2 > ${sample_id}_G_of_R_list.txt
  cp \$(ls ${sample_id}_L1_*.aln.bam_TE_list.txt | head -1) ${sample_id}_TE_list.txt #note: we only need one TE list file since they are all the same
  touch ${sample_id}_gofr.done
  """
}

process run_em {
  //maxForks 1  // limits to 1 concurrent execution
  tag "$sample_id"
  input:
    tuple val(sample_id), path(l1em_scripts), path(gofr_done), path(g_r_list), path(te_list), path(pk2s), path(tes), val(EM_threshold), val(nr_threads)

  output:
    tuple val(sample_id), path("${sample_id}_em.done"), path("${sample_id}_names_final.pkl"), path ("${sample_id}_X_final.pkl"), emit: final_x_pkl

  script:
  """
  python ${l1em_scripts}/L1EM/L1EM.py -g "${sample_id}_G_of_R_list.txt" -l ${sample_id}_TE_list.txt -t ${nr_threads} -s ${EM_threshold}
  mv names_final.pkl ${sample_id}_names_final.pkl
  mv X_final.pkl ${sample_id}_X_final.pkl
  touch "${sample_id}_em.done"
  """
}


// TO DO: revise all of these scripts for correctness and necessity. Right now, we only use the full_counts.txt to calculate L1 TPM values.
// TO DO: would be nice to add here a process to calculate the TPM values. For that, we need the RSEM quant outputs though for all the genes (or TPM values produced by other algorithm).
// TO DO: maybe use publishdir to directly put the results in the results_dir instead of copying them over at the end?

process finalize_results {
  tag "$sample_id"
  input:
    tuple val(sample_id), val(strandedness), path(l1em_scripts), val(results_dir), path(em_done), path(baminfo), path(g_r_list), path(pk2s), path(final_names_pkl), path(final_x_pkl)

  output:
    path("${sample_id}_l1em_done")

  script:
  """
  python ${l1em_scripts}/utilities/L1EM_readpairs_nf.py "${sample_id}_G_of_R_list.txt" >> ${sample_id}_baminfo.txt
  if [[ "$strandedness" == "stranded" ]]; then
    python ${l1em_scripts}/utilities/report_l1_exp_counts_nf.py "${sample_id}_G_of_R_list.txt" "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" > ${sample_id}_full_counts.txt
    python ${l1em_scripts}/utilities/report_l1hs_transcription_nf.py "${sample_id}_G_of_R_list.txt" "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" > ${sample_id}_l1hs_transcript_counts.txt
    python ${l1em_scripts}/utilities/filtered_and_normalized_l1hs.py "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" \$(head -2 ${sample_id}_baminfo.txt | tail -1) \$(head -3 ${sample_id}_baminfo.txt | tail -1) > ${sample_id}_filter_L1HS_FPM.txt
  else
    python ${l1em_scripts}/utilities/report_l1_exp_counts_unstranded_nf.py "${sample_id}_G_of_R_list.txt" "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" > ${sample_id}_full_counts.txt
    python ${l1em_scripts}/utilities/report_l1hs_transcription_unstranded_nf.py "${sample_id}_G_of_R_list.txt" "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" > ${sample_id}_l1hs_transcript_counts.txt
    python ${l1em_scripts}/utilities/filtered_and_normalized_l1hs_unstranded.py "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" \$(head -2 ${sample_id}_baminfo.txt | tail -1) \$(head -3 ${sample_id}_baminfo.txt | tail -1) > ${sample_id}_filter_L1HS_FPM.txt
  fi
    
  cp -rL ${sample_id}_baminfo.txt ${sample_id}_full_counts.txt ${sample_id}_l1hs_transcript_counts.txt ${sample_id}_filter_L1HS_FPM.txt "${sample_id}_names_final.pkl" "${sample_id}_X_final.pkl" ${results_dir}/
  touch "${results_dir}/${sample_id}_l1em_done"
  touch "${sample_id}_l1em_done" 
  """
}
