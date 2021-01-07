# Reactive astrocytes in ALS display dysregulated intron retention

This repository contains scripts to analyse the data and reproduce the figures from the paper:

[Reactive astrocytes in ALS display dysregulated intron retention](biorxiv)
Oliver J. Ziff, Doaa M. Taha, Hamish Crerar, Benjamin E. Clarke, Anob M. Chakrabarti, Gavin Kelly, Jacob Neeves, Giulia Tyzack, Nicholas M. Luscombe, Rickie Patani

The scripts are written in Rmarkdown documents for readability and are organised in order of the Figures in the paper.

All RNA sequencing data generated for this study is deposited at [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) under accession number [GSE160133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160133). RAW Mass Spectrometry data have been deposited to the ProteomeXchange Consortium (http://proteomecentral.proteomexchange.org) via the PRIDE partner repository with the dataset identifier PXD022604. 

Previously published iPSC derived astrocytes carrying ALS mutations are available at [GSE142730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142730) (C9orf72), [GSE102902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102902) and [GSE99843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99843) (SOD1 mutants and control respectively). Cytokine-stimulated iPSC derived astrocytes are available at [syn21861181](https://www.synapse.org/#!Synapse:syn21861181). TARDBP deleted mouse spinal cord astrocyte specific RNA-seq is available at [GSE156542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156542). Mouse SOD1 astrocyte TRAP-seq is available at [GSE74724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74724).

For each dataset we analyse:
1. QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks with [nf-core/rnaseq](https://github.com/nf-core/rnaseq) using the [MultiQC](https://multiqc.info/) html output
2. Alignment to genome with [STAR](https://github.com/alexdobin/STAR) and indexed with [SAMtools](https://sourceforge.net/projects/samtools/files/samtools/)
3. Read quantification using [HTSeq](https://htseq.readthedocs.io/en/master/)
4. Differential gene expression using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
5. Intron retention using [IRFinder](https://github.com/williamritchie/IRFinder)

STAR and samtools were run as follows:
```bash
STAR --runThreadN 1 --genomeDir $IDX --readFilesIn $READ1 $READ2 --outFileNamePrefix $OUT/${GROUP} --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --outSAMstrandField intronMotif
samtools index $OUT/${GROUP}Aligned.sortedByCoord.out.bam
```

HTSeq was run in intersection-strict mode as follows:
```bash
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $FILE $GTF > $OUT/${SAMPLE}.tab
```

DESeq2 was then run as per the DESeq2_analysis.Rmd script (see Figures folder).

IRFinder reference was built as follows:
```bash
IRFinder -m BuildRef -r $REF/Human-hg38-release99 ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
```

IRFinder was then run on merged technical repeat fastq files after adaptors had been trimmed using:

```bash
IRFinder -r $IRFINDER_REFERENCE -a none -d $OUT/$GROUP $READS
```

As we had two replicates we utilised the [audic and claverie test for differential IR](https://github.com/williamritchie/IRFinder/wiki/Small-Amounts-of-Replicates-via-Audic-and-Claverie-Test). We first pooled replicates of the same condition using IRFinder_pool_replicates scripts (see splicing folder). We then measured pooled IR in each condition using the IRFinder_differential script (see splicing folder).

Downstream analysis of DESeq2 and IRFinder outputs were then analysis as per the `figures_and_tables_resubmission.Rmd` script.

