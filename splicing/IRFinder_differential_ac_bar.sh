#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=IRFinComp
#SBATCH --output=IRFinComp-%j.out
#SBATCH --error=IRFinComp-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.ziff@crick.ac.uk

ml GCC Perl STAR BEDTools SAMtools
#ml GCC Perl STAR/2.7.1a-foss-2019a BEDTools SAMtools/1.9-foss-2018b  # load dependencies
OUT=/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder
cd $OUT

analysisWithLowReplicates.pl -A pooled_ac_bar_a1/IRFinder-IR-nondir.txt ac_bar_a1m/IRFinder-IR-nondir.txt ac_bar_a1f/IRFinder-IR-nondir.txt ac_bar_a1u/IRFinder-IR-nondir.txt -B pooled_ac_bar_a0/IRFinder-IR-nondir.txt ac_bar_a0m/IRFinder-IR-nondir.txt ac_bar_a0f/IRFinder-IR-nondir.txt ac_bar_a0u/IRFinder-IR-nondir.txt > ac_bar_a1-vs-a0.txt
