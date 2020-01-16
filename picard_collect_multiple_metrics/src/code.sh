#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail


#
# Fetch inputs
#
# ~/in/sorted_bam/*
#
mkdir -p ~/in/sorted_bam/
mkfifo "$sorted_bam_path"
dx cat "$sorted_bam" > "$sorted_bam_path" &
sorted_bam_pid=$!

#
# Stream and unpack genome
# 
# modified this from original to accept different genome reference file
#dx cat "$DX_PROJECT_CONTEXT_ID:/References/hs37d5.fa.gz" | gunzip -c > genome.fa
dx cat "$fasta_index" | tar zxvf - # => genome.fa, genome.fa.fai, genome.dict
#
# Set up options
#

# Calculate 80% of memory size, for java
mem_in_mb=`head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}'`
java="java -Xmx${mem_in_mb}m"

#
# Run Picard
#
mkdir -p out/stats/
$java -jar /picard.jar CollectMultipleMetrics I="$sorted_bam_path" R=genome.fa PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle PROGRAM=CollectBaseDistributionByCycle O="out/stats/$sorted_bam_prefix" $extra_options

wait "$sorted_bam_pid"

#
# Upload results
#
dx-upload-all-outputs --parallel
