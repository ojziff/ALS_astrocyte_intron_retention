#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=IRFinDiff
#SBATCH --output=IRFinDiff-%j.out
#SBATCH --error=IRFinDiff-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.ziff@crick.ac.uk

ml GCC Perl STAR/2.7.1a-foss-2019a BEDTools SAMtools/1.9-foss-2018b  # load dependencies
IRFINDER_REFERENCE=/camp/home/ziffo/home/genomes/ensembl/Human-hg38-release99
OUT=/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder

IRFinder -m BAM  -r $IRFINDER_REFERENCE -d $OUT/pooled_ac_cyt_ctrl <(samtools cat ac_cyt_ctrl1/Unsorted.bam ac_cyt_ctrl5/Unsorted.bam)
