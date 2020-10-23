#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=IRFinDiff
#SBATCH --output=IRFinDiff-%j.out
#SBATCH --error=IRFinDiff-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.ziff@crick.mn.uk

ml GCC Perl STAR BEDTools SAMtools  # load dependencies
IRFINDER_REFERENCE=/camp/home/ziffo/home/genomes/ensembl/Human-hg38-release99
OUT=/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder

IRFinder -m BAM  -r $IRFINDER_REFERENCE -d $OUT/pooled_ac_bar_a0 <(samtools cat ac_bar_a0m/Unsorted.bam ac_bar_a0f/Unsorted.bam ac_bar_a0u/Unsorted.bam)
