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

analysisWithNoReplicates.pl -A ac_tyz_sod1/IRFinder-IR-dir.txt -B ac_tyz_ctrl/IRFinder-IR-dir.txt > ac_tyz_sod1-vs-ctrl.txt
