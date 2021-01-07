#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=IRFinComp
#SBATCH --output=IRFinComp-%j.out
#SBATCH --error=IRFinComp-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.ziff@crick.ac.uk

ml GCC Perl STAR BEDTools SAMtools # load dependencies
OUT=/camp/home/ziffo/home/downloaded-public-data/astrocyte-tdp43-mouse-peng-2020/splicing/IRFinder
cd $OUT

analysisWithLowReplicates.pl -A pooled_ac_tdp_ko/IRFinder-IR-nondir.txt ac_tdp_ko1/IRFinder-IR-nondir.txt ac_tdp_ko2/IRFinder-IR-nondir.txt ac_tdp_ko3/IRFinder-IR-nondir.txt ac_tdp_ko4/IRFinder-IR-nondir.txt -B pooled_ac_tdp_ctrl/IRFinder-IR-nondir.txt ac_tdp_ctrl1/IRFinder-IR-nondir.txt ac_tdp_ctrl2/IRFinder-IR-nondir.txt ac_tdp_ctrl3/IRFinder-IR-nondir.txt ac_tdp_ctrl4/IRFinder-IR-nondir.txt > ac_tdp_ko_vs_ctrl.txt
