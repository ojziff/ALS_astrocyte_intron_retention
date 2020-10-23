#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=IRFinComp
#SBATCH --output=IRFinComp-%j.out
#SBATCH --error=IRFinComp-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.ziff@crick.ac.uk

ml GCC Perl STAR/2.7.1a-foss-2019a BEDTools SAMtools/1.9-foss-2018b  # load dependencies
OUT=/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder
cd $OUT

analysisWithLowReplicates.pl -A pooled_ac_cyt_vcp/IRFinder-IR-dir.txt ac_cyt_cb1e/IRFinder-IR-dir.txt ac_cyt_glia/IRFinder-IR-dir.txt -B pooled_ac_cyt_ctrl/IRFinder-IR-dir.txt ac_cyt_ctrl1/IRFinder-IR-dir.txt ac_cyt_ctrl5/IRFinder-IR-dir.txt > ac_cyt_vcp-vs-ctrl.txt
