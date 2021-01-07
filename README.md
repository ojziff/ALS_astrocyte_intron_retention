# Reactive astrocytes in ALS models display dysregulated intron retention

This repository contains scripts to analyse the data and reproduce the figures from the paper:

[Reactive astrocytes in ALS models display dysregulated intron retention](biorxiv)
Oliver J. Ziff, Doaa M. Taha, Hamish Crerar, Benjamin E. Clarke, Anob M. Chakrabarti, Gavin Kelly, Jacob Neeves, Giulia Tyzack, Nicholas M. Luscombe, Rickie Patani

The scripts are written in Rmarkdown documents for readability and are organised according to the Figures in the paper.

Raw data are deposited at [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) under accession number [GSE160133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160133).

Previously published iPSC derived astrocytes carrying ALS mutations are available at [GSE142730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142730) (C9orf72), [GSE102902](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102902) and [GSE99843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99843) (SOD1 mutants and control respectively). Cytokine-stimulated iPSC derived astrocytes are available at [syn21861181](https://www.synapse.org/#!Synapse:syn21861181). TARDBP knockout mouse spinal cord astrocyte specific RNA-seq is available at [GSE156542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156542). Mouse SOD1 astrocyte TRAP-seq is available at [GSE74724](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74724).

For each dataset we analyse:
1. Differential gene expression using HTSeq and DESeq2
2. Intron retention using IRFinder

HTSeq was run in intersection-strict mode as follows:
```bash
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $FILE $GTF > $OUT/${SAMPLE}.tab
```

DESeq2 was then run as per the DESeq2_analysis.Rmd script (see gene_expression folder).

IRFinder reference was built as follows:
```bash
IRFinder -m BuildRef -r $REF/Human-hg38-release99 ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz
```

IRFinder was then run on merged technical repeat fastq files after adaptors had been trimmed using:

```bash
IRFinder -r $IRFINDER_REFERENCE -a none -d $OUT/$GROUP $READS
```

As we had two replicates we utilised the [audic and claverie test for differential IR](https://github.com/williamritchie/IRFinder/wiki/Small-Amounts-of-Replicates-via-Audic-and-Claverie-Test). We first pooled replicates of the same condition using IRFinder_pool_replicates scripts (see splicing folder). We then measured pooled IR in each condition using the IRFinder_differential script (see splicing folder).

Downstream analysis of DESeq2 and IRFinder outputs were then analysis as per the figures_and_tables script.

