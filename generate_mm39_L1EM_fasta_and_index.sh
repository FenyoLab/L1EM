#! /bin/bash

# If you need to specify package directories
bedtools=$(which bedtools)
bwa=$(which bwa)

# Command line
mm39=$1

$bedtools getfasta -s -name -fi $mm39 -bed annotation/mm39.L1EM.bed > annotation/mm39.L1EM.400.fa
$bwa index annotation/mm39.L1EM.400.fa
