#!/bin/bash -l

#SBATCH -t 48:00:00
#SBATCH --mem=40G
#SBATCH -c 10
#SBATCH --job-name=IRFinDiff
#SBATCH --output=IRFinDiff-%j.out
#SBATCH --error=IRFinDiff-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=oliver.ziff@crick.mn.uk

ml GCC Perl STAR BEDTools SAMtools/1.9-foss-2018b  # load dependencies
IRFINDER_REFERENCE=/camp/home/ziffo/home/genomes/ensembl/Mouse-mm38-release101
OUT=/camp/home/ziffo/home/downloaded-public-data/astrocyte-tdp43-mouse-peng-2020/splicing/IRFinder

IRFinder -m BAM  -r $IRFINDER_REFERENCE -d $OUT/pooled_ac_tdp_ctrl <(samtools cat ac_tdp_ctrl1/Unsorted.bam ac_tdp_ctrl2/Unsorted.bam ac_tdp_ctrl3/Unsorted.bam ac_tdp_ctrl4/Unsorted.bam)
