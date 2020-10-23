.libPaths("/camp/lab/luscomben/home/.conda/envs/rtest/lib/R/library")
# library(devtools)
# install_github("jsha129/irlib")
library(irlib)
library(tidyverse)
library(DESeq2)
library(BiocParallel)
library("grid")
library("ggplotify")
library("RColorBrewer")
library("ggpubr")
library("gprofiler2")
library(readxl)
library(rhdf5)
library(dplyr)
library(purrr)
library("stringr")
library("tximeta")
library("tximport")
library(rtracklayer)
library("DESeq2")
library("rmarkdown")
library("data.table")
library("gsubfn")
library('ggplot2')
library("ggpubr")
library("pheatmap")
library('ggrepel')
library("grid")
library("ggplotify")
library("clusterProfiler")
library("BiocParallel")
library("EnhancedVolcano")
library("org.Hs.eg.db")
library("limma") 
library("AnnotationDbi")
library("Glimma")
library("colorRamps")
library("RColorBrewer")
library("GO.db")
library("fgsea")
library("geneplotter")
library("genefilter") 
library("sva")
library("systemPipeR")
library("vsn")
library("RUVSeq")
library("gprofiler2")
library("topGO") 
library(Biobase)
library(SummarizedExperiment)
library(VennDiagram)
library(ggplotify)
library(patchwork)
library(rstatix)
library(data.table)
library(GenomicFeatures)
library(Gviz)
library(trackViewer)
library(GeneOverlap)
library(ggdendro)
library(dendextend)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Rnorvegicus.UCSC.rn6)
library(GenomicScores)
library(phastCons100way.UCSC.hg38)
library(DEP)
library(ggforce)
library(pcr)
library(gplots)
library(notNMD)
library(plyranges)
library(seqinr)
library(GeneStructureTools)
library(GenomicRanges)
library(notNMD) # devtools::install_github('betsig/notNMD')
library(gdata)
library(DEP)
options(stringsAsFactors = F)

# GO terms ----------------------
gprofiler_database <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/gprofiler_full_hsapiens.name.edit.gmt") # IDs removed to see term names
# https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp GO pathways
curated.gene.sets.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c2.all.v7.1.symbols.gmt") # includes KEGG, REACTOME, PID, BIOCARTA, CGP
kegg.pathways.msigdb<- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c2.cp.kegg.v7.1.symbols.gmt")
reactome.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c2.cp.reactome.v7.1.symbols.gmt")
go.pathways.msigdb<- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.all.v7.1.symbols.gmt")
go.bp.pathways.msigdb<- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.bp.v7.1.symbols.gmt")
go.mf.pathways.msigdb<- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.mf.v7.1.symbols.gmt")
go.cc.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c5.cc.v7.1.symbols.gmt")
immune.pathways.msigdb <- gmtPathways("/camp/home/ziffo/home/genomes/gene-ontology/c7.all.v7.1.symbols.gmt")
RNA_BINDING_PROTEINS <- read_excel("/camp/home/ziffo/home/genomes/gene-ontology/rbp_consensus_gerstberger_2014.xls", sheet = "RBP table") %>% pull(`gene name`)
glia.genes <- c("S100B", "SOX9")
astrocyte.genes <- c("GFAP","AQP4","PLA2G7","SLC39A12","MLC1","DIO2","SLC14A1","ALDH1L1","ALDOC","TTPA","ACSBG1","CHRDL1","SLC4A4","SLC1A2","SLC25A18","SLC1A3","F3","PPP1R3G","FZD2","MERTK","EZR","EVA1A","GJB6","HAPLN1","RFX4","PAPSS2","SLC15A2","PPP1R3C","TLR3","ACOT11","ATP1A2","BMPR1B","PRODH","GLI3","TMEM47","SLC9A3R1","CTH","NTSR2","SLC7A10","VCAM1","FGFR3","CCDC80","ENTPD2","CYBRD1","KCNE5","FAM20A","TNC","TLCD1","S1PR1","CBS","PBXIP1","GRIN2C","ADHFE1","AGT","GLDC","SLC7A2","GJA1","PDK4","EGFR","SOX9","CLDN10","PLCD4","ID4","FMO1","EMP2","LONRF3","HTRA1","MGST1","THRSP")
panreactive.astrocyte.genes <- c("LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "OSMR", "CP", "VIM", "GFAP",
                                 "SERPINA3", "ASPG") # "SERPINA5",
A1.astrocyte.genes <- c("C3", "HLA-A", "HLA-B", "HLA-C","HLA-E","HLA-F", "MICA", "H2-T23","H2-D1", #, "HCP5", "HLA-H", "HLA-G", "HLA-K", "HLA-L", "AL645929.2",
                        "GBP2", "AMIGO2", "SERPING1","GGTA1P","GGTA1", #"GLT6D1", "A3GALT2", "IRGC",
                        "FBLN5", "UGT1A1", "FKBP5", "PSMB8", "SRGN","IIGP1")#, "MX1")
A2.astrocyte.genes <- c("S100A10", "EMP1", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")
astrocyte.subtype.genes <- list(panreactive.astrocyte.genes,A1.astrocyte.genes,A2.astrocyte.genes)
astrocyte.reactive.genes <- c(panreactive.astrocyte.genes, A1.astrocyte.genes, A2.astrocyte.genes) %>% unique
astrocyte.reactivity.markers = data.frame("gene_name" = astrocyte.reactive.genes) %>% 
  mutate(group = case_when(gene_name %in% A1.astrocyte.genes ~ "A1",
                           gene_name %in% A2.astrocyte.genes ~ "A2",
                           gene_name %in% panreactive.astrocyte.genes ~ "Pan-reactive"))
all.astrocyte.markers <- c("GRN", "PSEN1", "LRP1", "APP", "VIM", "IFNGR1", "GRN", "LAMC3", "PSEN1", "MT3", "TREM2", "CDK6", "IFNG", "LRP1", "PLP1", "IL1B", "NR1D1", "ADORA2A", "SMO", "LDLR", "GFAP", "KRAS", "TSPAN2", "IL6", "TLR4", "APP", "S100A8", "EIF2B5", "TTBK1", "EGFR", "ROR2", "FPR2", "LAMB2", "C1QA", "BACE2", "POU3F2", "DRD1", "ROR1", "MAPT", "NF1", "C5AR1", "DLL1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "KRAS", "NF1", "CCR2", "GPR183", "WNT1", "FZD1", "CTNNB1", "HEXB", "CCL2", "CCR2", "APCDD1", "MMP14", "GPR183", "SCRIB", "CCL3", "IFNGR1", "GRN", "PSEN1", "TREM2", "IFNG", "LRP1", "IL1B", "NR1D1", "ADORA2A", "SMO", "LDLR", "IL6", "APP", "TTBK1", "EGFR", "FPR2", "C1QA", "BACE2", "MAPT", "C5AR1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "SOX8", "PAX6", "VIM", "IFNGR1", "GRN", "LAMC3", "PSEN1", "MT3", "TREM2", "ABL1", "MAPK1", "MAPK3", "MAG", "CDK6", "PRPF19", "SOX6", "IFNG", "NR2E1", "HES1", "ID2", "EPHA4", "LRP1", "PLP1", "TTC21B", "SOX9", "IL1B", "NKX2-2", "BMP2", "NR1D1", "ADORA2A", "LIF", "SMO", "LDLR", "GFAP", "KRAS", "TSPAN2", "IL6ST", "SERPINE2", "IL6", "BIN1", "TLR4", "GCM1", "NTRK3", "MBD1", "APP", "S100A8", "EIF2B5", "TTBK1", "EGFR", "NOTCH1", "S100B", "TAL1", "SHH", "STAT3", "MAP2K1", "ROR2", "GPR37L1", "FPR2", "LAMB2", "ID4", "C1QA", "DAB1", "CLCF1", "PTPN11", "F2", "BACE2", "NOG", "CNTN2", "POU3F2", "DRD1", "ROR1", "MAPT", "NF1", "C5AR1", "HES5", "DLL1", "AGER", "MIR181C", "MIR181B2", "MIR181B1", "TNF", "CNTF", "MIR142", "MAG", "PRPF19", "NR2E1", "HES1", "ID2", "EPHA4", "BMP2", "NR1D1", "LIF", "LDLR", "IL6ST", "SERPINE2", "IL6", "BIN1", "NTRK3", "MBD1", "TTBK1", "NOTCH1", "GPR37L1", "ID4", "DAB1", "CLCF1", "F2", "NOG", "CNTN2", "NF1", "HES5", "MIR181C", "MIR181B2", "MIR181B1", "MIR142", "MAG", "PRPF19", "HES1", "ID2", "BMP2", "LIF", "IL6ST", "SERPINE2", "BIN1", "TTBK1", "NOTCH1", "CLCF1", "MIR142", "NR2E1", "NR1D1", "LDLR", "NTRK3", "MBD1", "GPR37L1", "ID4", "DAB1", "F2", "NOG", "NF1", "HES5", "MIR181C", "MIR181B2", "MIR181B1", "SOX8", "SOX9", "GCM1", "TAL1", "NR1D1", "LDLR", "IL6", "TTBK1", "MIR181C", "MIR181B2", "MIR181B1", "MIR142", "NR1D1", "LDLR", "MIR181C", "MIR181B2", "MIR181B1", "TTBK1", "MIR142", "CCR2", "GPR183", "KCNK2", "MT3", "EZR", "MLC1", "SLC1A2", "ATP1B2", "GFAP", "SYT4", "EIF2S1", "APP", "SLC7A11", "PINK1", "GRM2", "GJB2", "AQP4", "KCNJ10", "SLC17A8", "GRM3", "DMD", "ADGRG1", "MT3", "MLC1", "ATP1B2", "GFAP", "EIF2S1", "AQP4", "SLC17A8", "ADGRG1", "MLH1", "MSH2", "TSC2", "IFNG", "MSH3", "MSH6", "PMS2", "POT1", "APC", "IDH1", "BRCA2", "TP53", "ERBB2", "CDKN2A", "AIFM1", "TSC1", "IDH2", "NF2", "NF1")
all.astrocyte.markers <- c(all.astrocyte.markers,astrocyte.reactive.genes)
microglia.genes <- c("ITGAM", "CX3CR1", "CCL3", "CSFIR", "CCL4", "P2RY12", "C1QB", "PLEK", "GPR183")
focal.adhesion.go <- c("LAP3", "CD99", "LASP1", "ABCB4", "ITGA3", "ITGA2B", "RALA", "CD9", "MRC2", "TSPAN9", "PLAUR", "EHD3", "CAPN1", "FHL1", "VIM", "CD44", "ARHGAP31", "VCL", "TNC", "CTNNA1", "HSPA5", "LIMA1", "BCAR1", "CYBA", "SYNE2", "GDI2", "PPP1R12A", "NCKAP1", "RPL18", "CNN2", "SLC9A3R2", "TLE2", "RHOA", "FGFR3", "PABPC1", "CDC42", "MAP4K4", "RPL31", "ACTN1", "LIMS2", "PVR", "FERMT2", "CLASP1", "HACD3", "ACTB", "REXO2", "MCAM", "USP33", "APBB1IP", "ACTN2", "ITGA8", "FAP", "TNS1", "SENP1", "DNM2", "EPB41L2", "RAB21", "PTPRC", "ITGB5", "RPS5", "FAT1", "RAB10", "CD59", "CPNE3", "CTTN", "NOX4", "TRIP6", "ADD1", "CASS4", "ASAP3", "RPL6", "RPLP0", "PXN", "SLC9A1", "ICAM1", "ITGA6", "SNAP23", "EZR", "SORBS1", "JAK2", "NRP1", "MISP", "MAPK1", "TRIOBP", "PACSIN2", "RPL3", "MYH9", "ZFYVE21", "PROCR", "FERMT1", "HCK", "MAPRE1", "CD99L2", "ARHGEF7", "FLT1", "MAPK3", "CORO2B", "RPS16", "RPS19", "ITGB8", "CAV2", "CAV1", "HSPB1", "LIMK1", "ENG", "PDLIM1", "GIT1", "RPL19", "PFN1", "SLC6A4", "YWHAE", "LAMTOR3", "HSPA8", "LPXN", "CBL", "CD81", "RPS13", "PPFIBP1", "CORO1C", "TNS2", "TRPV4", "ARPC3", "NEDD9", "OPRM1", "WASF1", "PTK7", "HSPA9", "PDGFRB", "PRKAR2A", "ACTR3", "EPB41L5", "ITGB6", "ITGA4", "RPS15", "IL1RL1", "FHL2", "RND3", "RPL22", "RHOU", "ARHGEF2", "DOCK7", "CD46", "CNN3", "GNA13", "TEK", "SORBS3", "PTK2B", "CAT", "RPL5", "PLAU", "ADGRE5", "LRP1", "SDC4", "STX16", "RPS10", "AHNAK", "EFNB2", "RPL23", "VASP", "FLRT3", "RRAS", "AP006333.1", "AIF1L", "MAP2K2", "PTPN12", "CDC42EP1", "RAC2", "FLNC", "ARHGAP22", "PALLD", "AJUBA", "STARD8", "CNN1", "NECTIN2", "ACTN4", "ARPC1B", "PAK4", "AKAP12", "CAP1", "RPL27", "PPFIA1", "TNS4", "ITGB4", "FLOT2", "PTPRA", "DCTN4", "MPRIP", "KRAS", "RRAS2", "NUMB", "YWHAQ", "ANXA1", "TES", "AVIL", "FLNB", "LMO7", "LCP1", "TNS3", "RAC1", "ARPC5L", "TLN1", "HMGA1", "FLOT1", "MDC1", "ATAT1", "SDCBP", "RDX", "KIF23", "ITGA11", "RPLP1", "ADAM10", "BCAR3", "ACTR2", "ITGAV", "ARHGAP24", "SCARB2", "PARVG", "GIT2", "ITGB7", "IQGAP1", "TGFB1I1", "ARMC5", "CDH13", "RPS2", "PDPK1", "CLTC", "GRB7", "RPS11", "RPL13A", "EPHA2", "HSPG2", "RPS8", "MTF2", "DCAF6", "PRUNE1", "PIP5K1A", "S100A7", "ARF1", "RHOB", "FBLN7", "LIMD1", "PHLDB2", "LPP", "RPS3A", "ARHGAP26", "G3BP1", "GNA12", "EGFR", "SH3KBP1", "CASK", "MSN", "ZNF185", "RPL7", "GSN", "RPL7A", "RSU1", "CAPN5", "PAK1", "RPS3", "HYOU1", "ITGB1", "DIXDC1", "TWF1", "ADAM17", "DST", "ARL14EP", "TADA1", "GJA1", "DAB2", "THY1", "PGM5", "ENAH", "SORBS2", "RPL30", "MMP14", "FZD1", "SHROOM4", "CSRP1", "ACTC1", "ZYX", "ITGB2", "RPL8", "NPHS1", "ITGA5", "JAK1", "FBLIM1", "NEXN", "ARPC5", "DDR2", "NCSTN", "CAPN2", "XIRP2", "ARPC2", "NFASC", "CLASP2", "RPL9", "ANXA5", "ITGA2", "RPS14", "DLC1", "SLC4A2", "YWHAZ", "HNRNPK", "ARF6", "ILK", "HSP90B1", "B2M", "PPIB", "YWHAB", "MAPRE2", "PDIA3", "TPM4", "SRP68", "CTNNB1", "FAM107A", "IRF2", "XIRP1", "ADAM9", "SNTB2", "TM4SF20", "MAP2K1", "PTK2", "LIMS1", "ALCAM", "YWHAG", "PDCD6IP", "CDH2", "RPS9", "RPS7", "TLN2", "KLF11", "SNTB1", "BSG", "GNB2", "SYNPO2", "CORO1B", "CFL1", "RPL38", "DAG1", "PEAK1", "CSPG4", "JUP", "RPL4", "CSRP2", "YES1", "RHOG", "RPLP2", "CD151", "FLII", "PLEC", "GAK", "CALR", "FZD2", "NPM1", "ADGRB1", "FES", "CAV3", "RPS17", "FHL3", "ACTG1", "UBOX5", "FLRT2", "LMLN", "P4HB", "ATP6V0C", "PIP5K1C", "PPP1CC", "CHP1", "SPRY4", "NHS", "PEAK3", "FOCAD", "PARVB", "PPIA", "EVL", "AFAP1", "MME", "PDLIM7", "FLNA", "ANXA6", "IGF2R", "PCBP2", "SRC", "SVIL", "DPP4", "PARVA", "RPL37A", "RPL12", "MPZL1", "RPS4X", "ITGBL1", "RPL10A", "L1CAM", "TGM2", "LAYN", "HSPA1B", "HSPA1A", "ARL2", "PPP1CB", "RPS29", "ITGA1", "TSPAN4", "RPS18", "ALKBH6", "PI4KA", "SCARF2", "ACTN3", "LIMS4", "LIMS3", "ITGB3", "AC068234.1", "CYFIP1", "PRAG1", "MARCKS")
collagen.formation.go <- c("FOXC1","COL11A1","COL5A3","ADAMTS2","TGFB2","MMP11","CHADL","COMP","PLOD3","AEBP1","TGFBR1","COL1A1","COL12A1","LOX","LOXL3","FMOD","P4HA1","COLGALT1","PXDN","COL5A1","LOXL2","CYP1B1","EMILIN1","LOXL4","ADAMTS14","COL2A1","LUM","RB1","P3H4","DPT","SFRP2","SERPINH1","VIPAS39","ADAMTS3","ACAN","DDR2","COL1A2","ATP7A","GREM1","SERPINF2","TNXB","COL3A1","FOXC2","ANXA2","COL14A1","OPTC","NF1","COL11A2","COL5A2","SCX","MIR29B1","MIR29B2 MMP25","MRC2","FAP","COL19A1","ADAMTS2","MMP2","MMP11","MMP9","CST3","MMP15","VSIR","CTSD","MMP8","MMP19","PEPD","MMP24","CTSL","MMP7","MMP20","MMP27","MMP13","ADAMTS14","FURIN","CTSK","ADAM15","MMP3","ITGB1","MMP21","MMP16","ADAMTS3","MMP14","CTSS","CTSB","MMP10","MMP26","KLK6","PHYKPL","TMPRSS6","MMP23B","PRTN3","MMP1","COL13A1","MMP17","COL15A1","MMP12","MMP28","PRSS2","MMP25","MRC2","VIM","TRAM2","FAP","COL19A1","ADAMTS2","MMP2","P2RX7","P3H2","SUCO","MMP11","HIF1A","MMP9","CST3","RGCC","MMP15","SMPD3","TGFB1","PLOD3","ENG","VSIR","COL1A1","P3H3","TNS2","PPARD","PDGFRB","ERRFI1","RAP1A","P3H1","CTSD","MMP8","MYB","ARG1","CCN2","TGFB3","GOT1","MMP19","PEPD","AMELX","BMP4","MMP24","ID1","COL5A1","PPARG","CTSL","IL6","MMP7","MMP20","MMP27","MMP13","EMILIN1","ADAMTS14","FURIN","ARRB2","CBX8","P3H4","RCN3","CTSK","ADAM15","SERPINH1","MMP3","ITGB1","VIPAS39","MMP21","MMP16","ADAMTS3","MMP14","CREB3L1","CYGB","WNT4","CTSS","NPPC","IHH","UCN","ITGA2","COL1A2","CTSB","LARP6","SERPINB7","MFAP4","MMP10","MMP26","SERPINF2","KLK6","TNXB","PHYKPL","CIITA","F2","F2R","VPS33B","TMPRSS6","MMP23B","PRTN3","HDAC2","MMP1","COL13A1","MMP17","COL15A1","MIR149","MIR218-1","MIR218-2","SCX","MMP12","MMP28","PRSS2","MIR145","MIR92A1","MIR29B1","MIR29A","MIR21","MIR29B2","MIR92A2","COL4A4","COL1A1","COL4A2","ITGA11","UBASH3B","DDR2","ITGA2","SYK","COL4A3","OSCAR","COL4A1","COL4A5","COL4A6","DDR1","TLL1","COL9A2","COL23A1","LAMA3","LAMC2","COL11A1","COL17A1","P4HA2","COL5A3","COL4A4","COL19A1","PLOD1","COL16A1","ADAMTS2","P3H2","ITGA6","COL9A3","TLL2","MMP9","COL20A1","PCOLCE","PLOD3","COL1A1","P3H3","COL12A1","COL9A1","LOX","COL7A1","LOXL3","P3H1","P4HA1","COL10A1","COL21A1","LOXL1","COLGALT1","PXDN","COL5A1","ITGB4","LOXL2","COL4A2","CTSL","CTSV","MMP7","MMP20","MMP13","LOXL4","ADAMTS14","COL2A1","COL6A1","COL6A2","COL8A1","SERPINH1","P4HA3","MMP3","PLOD2","ADAMTS3","COL26A1","CTSS","COL6A3","PCOLCE2","COL1A2","CTSB","PPIB","BMP1","COL3A1","COL4A3","COL22A1","CRTAP","COL24A1","COL8A2","COL6A5","CD151","PLEC","COL18A1","P4HB","COL4A1","COL14A1","COL4A5","COL25A1","COL27A1","LAMB3","COL13A1","COL4A6","COLGALT2","COL11A2","COL5A2","COL15A1","COL6A6","COL28A1")
extracellular.matrix.go <- c("CFLAR","ST7","ITGAL","ITGA3","ITGA2B","NOX1","MMP25","DCN","CAPN1","PHLDB1","CD44","TNFRSF1B","IBSP","TIMP2","TLL1","VCAN","CDH1","TNC","GPM6B","COL9A2","ELN","LAMC3","COL23A1","LAMA3","FOXC1","LAMC2","COL11A1","COL17A1","TNFRSF1A","BCL3","CLASP1","NTN4","FSCN1","ICAM3","NFKB2","FBLN1","ITGA8","FAP","COL5A3","CPB2","COL4A4","COL19A1","ITGB5","ITGAE","COL16A1","B4GALT1","ERO1B","ADAMTS2","MMP2","NID2","KIF9","ICAM1","LAMB4","LAMB1","ITGA6","CCDC80","CMA1","COL9A3","TGFB2","TMEM38B","TLL2","ABL1","MADCAM1","MMP11","PDGFB","CHADL","CTSG","MMP9","FERMT1","CST3","LAMA1","TIMP1","RGCC","MMP15","HAS3","SMPD3","CRISPLD2","FOXF1","ERCC2","TGFB1","ICAM4","ICAM5","HAS1","COMP","HPN","ITGB8","CAV2","CAV1","DNAJB6","SERPINE1","PLOD3","AEBP1","MEGF9","TGFBR1","ECM2","ENG","SPOCK2","KAZALD1","SH3PXD2A","ICAM2","COL1A1","VTN","VWF","MYF5","COL12A1","ADTRP","COL9A1","NR2E1","SMOC2","LAMA4","LOX","SPARC","COL7A1","MPV17","ITGB6","ITGA4","LOXL3","FN1","TNR","QSOX1","EXOC8","NID1","MFAP2","MMP8","TTR","CCN2","SPP1","TGFBI","FMOD","PLG","SERAC1","RECK","P4HA1","PRDX4","MMP19","LRP1","COL10A1","MATN4","SOX9","TCF15","MMP24","CAPNS1","KDR","LOXL1","NCAN","COLGALT1","PXDN","COL5A1","LAMA5","POMT1","GFAP","RAMP2","MATN3","ITGB4","BCAN","HAPLN2","POSTN","MYH11","SPINK5","LOXL2","PDGFRA","COL4A2","ETS1","CTSL","ADAM19","ITGA7","AGT","LAMC1","LCP1","IL6","CTSV","FOXF2","FLOT1","SULF1","MMP7","MMP20","MMP27","MMP13","THBS1","ITGA11","ADAM10","CYP1B1","EMILIN1","LOXL4","ADAMTS14","ITGAV","FGF2","SLC39A8","FBN2","COL2A1","LUM","ITGB7","RB1","FBLN5","FURIN","ITGAX","GFOD2","P3H4","COL6A1","COL6A2","APP","HSPG2","CCN1","ITGA10","DPT","ADAMTSL4","CTSK","ADAM15","ITGA9","COL8A1","PHLDB2","SFRP2","HAPLN1","SCUBE3","CSGALNACT1","NOTCH1","ADAM12","GAS2","HSD17B12","SERPINH1","MMP3","ITGB1","VIPAS39","ADAM8","DSPP","DMP1","ABI3BP","WNT3A","MMP21","JAM2","ADAMTS5","MMP16","ADAMTS3","ITGAD","MMP14","MYO1E","CREB3L1","ACAN","F11R","ADAMTS4","SCUBE1","CARMIL2","ITGB2","MPZL3","FGFR4","NPHS1","ITGA5","PDPN","MATN1","VCAM1","DDR2","OLFML2B","CAPN2","CTSS","COL6A3","ELF3","IHH","FBLN2","CLASP2","ADAMTS9","PTX3","APBB2","MELTF","ITGA2","EGFLAM","KLKB1","COL1A2","TNFRSF11B","ATP7A","HTRA1","JAM3","SPINT1","FBN1","WDR72","MFAP4","MMP10","GREM1","SMAD3","MMP26","SPINT2","SERPINF2","KLK4","KLK2","KLK5","TNXB","BMP1","COL3A1","NPNT","CTRB1","CTRB2","COL4A3","KLK7","PTK2","COL22A1","ANTXR1","ITGAM","FSHR","HAS2","COL24A1","RXFP1","FGG","FGA","FGB","COL8A2","ANGPTL7","LAMB2","TPSAB1","BSG","EFEMP2","HPSE2","ADAMTS20","NDNF","DAG1","SH3PXD2B","A2M","FOXC2","RIC8A","VWA1","WASHC1","OTOL1","BGN","ANXA2","COL18A1","GAS6","WT1","FLRT2","OLFML2A","TMPRSS6","COL4A1","THSD4","COL14A1","COL4A5","AGRN","OPTC","MMP23B","NOXO1","SULF2","LAMA2","MMP1","NF1","COL27A1","CD47","LAMB3","PDGFA","COL13A1","ELANE","COL4A6","MFAP5","DPP4","ADAMTSL2","ERO1A","MMP17","EGFL6","COL11A2","COL5A2","COL15A1","DDR1","PRSS1","VIT","SERPINB5","COLQ","ITGA1","COL28A1","ATXN1L","TNF","MARCOL","CAPNS2","ITGB3","SCX","PECAM1","MMP12","MMP28","MIR98","PRSS2","MIR29B1","MIR29B2")
all.reactivity.go <- c(focal.adhesion.go, collagen.formation.go, extracellular.matrix.go, go.pathways.msigdb$GO_IMMUNE_EFFECTOR_PROCESS, gprofiler_database$`Cellular responses to stress`, all.astrocyte.markers)
astrocyte.reactivity.splicing.labels = c("VIM", "COL1A1", "MYL6", "DDX5", "FN1", "FLNA", "FUS", "SFPQ","OSMR", "NONO", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3","UBN1", "COL7A1","MX1",
                                "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4", "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2", "TTBK1", "CCS", "TGFB1I1","COL1A1A", "TNC", 
                                "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1")
astrocyte.reactivity.labels = c("VIM", "COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3","UBN1", "COL7A1","MX1",
                              "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4", "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2", "TTBK1", "CCS", "TGFB1I1","COL1A1A", "TNC", 
                              "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1",
                              "TINF2", "ILK", "MCAM", "FN1", "FBLN5", "CAPN2", "TGFB1I1", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1", "BAG3", "PVR")  # overlap increased exp & decreased IR
astrocyte.reactivity.labels = c(astrocyte.reactivity.labels, astrocyte.reactive.genes)

astrocyte.reactivity.labels.ac_nuc <- c("VIM", "COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6",
                                         "PIMREG", "RRP8", "RECK", "MX1", "TTBK1", "MBD1", "EMP1", "PRPF4", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ANKZF1", "GTBP2",   
                                         "CCS", "TGFB1I1", "COL1A1A", "RPL10", "NUP199", "LAMB2", "IDH1", "LRP1", "TINF2", "ILK" , "MCAM", "CAPN2", "ADAM19", "SERPINE1", "ITGA6", "PDLIM7", "FLOT1",  
                                         "BAG3", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3", "ASPG", "C3", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  
                                         "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_cyt <- c("VIM", "COL1A1", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "LOXL3", "COL7A1", "MX1", "EMP1",
  "NARPT", "CSRP1", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2",   "CCS", "TGFB1I1", "COL1A1A", "TNC", "NUP199", "LAMB2", "IDH1", "TP53", "ILK" ,    
  "MCAM", "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1",  "BAG3", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_nucIR_cytDGE <- c("VIM", "COL1A1", "FN1", "HLA-A", "HLA-B", "HLA-C", "SLC26AA6","PIMREG", "RRP8", "RECK", "MX1", "TTBK1", "EMP1", "PRPF4",   
  "NARPT", "CSRP1", "IRF7", "GBP2", "ANKZF1", "GTBP2","CCS", "TGFB1I1", "COL1A1A", "RPL10", "NUP199", "TP53", "TINF2", "ILK",    
  "MCAM", "CAPN2", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "FLOT1","BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_c9orf72 <- c("VIM", "COL1A1", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6","PIMREG", "RRP8", "RECK", "UBN1", "MX1", "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4",   
  "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2","CCS", "TGFB1I1", "COL1A1A", "TNC", "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1", "TINF2", "ILK" ,    
  "MCAM", "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "FLOT1", "DIXDC1", "BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "HLA-F", "MICA", "H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "SPHK1", "CD109", "PTGS2", "SLC10A6", 
  "TM4SF1", "B3GNT5", "CD14")  

astrocyte.reactivity.labels.ac_sod1 <-  c("COL1A1", "MYL6", "FN1", "FLNA", "OSMR", "HLA-A", "HLA-B", "HLA-C", "STAT3", "RND3", "SLC26AA6",
  "PIMREG", "RRP8", "RECK", "LOXL3", "UBN1", "COL7A1", "MX1", "TTBK1", "TSC1", "MBD1", "EMP1", "PRPF4",   
  "PIF1", "AARS2", "PLSCR1", "NARPT", "CSRP1", "IRF7", "ID4", "AGER", "GBP2", "ANKZF1", "AGRN", "GTBP2",   
  "CCS", "TGFB1I1", "COL1A1A", "TNC", "RPL10", "NUP199", "LAMB2", "IDH1", "TP53", "LRP1", "TINF2", "ILK" ,    
  "FBLN5", "CAPN2", "CDKN1A", "ADAM19", "SERPINE1", "ITGA6", "ADAMTSL4", "ITGA5", "PDLIM7", "DIXDC1",  
  "BAG3", "PVR", "LCN2", "STEAP4", "S1PR3", "TIMP1", "HSPB1", "CXCL10", "CD44", "CP", "GFAP", "SERPINA3",
  "ASPG", "C3", "HLA-E", "HLA-F","H2-T23", "H2-D1", "AMIGO2", "SERPING1", "GGTA1P", "GGTA1", "UGT1A1",  
  "FKBP5", "PSMB8", "SRGN", "IIGP1", "S100A10", "CLCF1", "TGM1", "PTX3", "SPHK1", "CD109", "PTGS2", "SLC10A6", 
  "TM4SF1", "B3GNT5", "CD14")  
astrocyte.reactivity.labels.ac_a1 <- c("FN1", "ILK", "OSMR", "TINF2", "IRF7", "MCAM", "SOD2", "PSMB7", "PSMB10", "DNAJC2", "KMT2E", "COL27A1", "ERF")

RNA.capping <-c("POLR2J","MNAT1","POLR2B","NCBP3","POLR2E","POLR2F","RNMT","POLR2C","ERCC2","POLR2I","GTF2H1","GTF2H3","RNGTT","NCBP2","GTF2F1","CDK7","CCNH","NCBP1","CMTR1","TGS1","POLR2D","GTF2H2","POLR2K","ERCC3","POLR2H","POLR2G","RAMAC","POLR2L","CMTR2","POLR2A","GTF2F2","SUPT5H","GTF2H4","RAMACL","GTF2H5")
exon.junction.complex <- c("UPF1","THRAP3","SMG6","TDRD3","PNN","R3HCC1","CASC3","MAGOHB","UPF3B","SRSF1","EIF4A3","SAP18","UPF2","MAGOH","R3HCC1L","UPF3A","PYM1","RNPS1","RBM8A","SRRM1","PRPF8","EIF4A3","MLN51","MAGOH","Y14","CASC3","EIF4A3","MAGOH","RBM8A","PYM1","NXF1","DDX39B")
RNA.polyadenylation <- c("PAF1","ZC3H3","PABPC1","CPSF1","SNRPA","YTHDC1","PAPOLA","ZC3H14","PABPC1L","CSTF2","RNF40","CPEB3","MTPAP","CPSF6","TENT4A","PAPOLG","CPSF3","TENT4B","SYMPK","CCNT1","GRSF1","CDC73","WDR33","CDK9","PNPT1","APP","FIP1L1","TUT1","CPSF7","RNF20","SSU72","TENT2","VIRMA","PCF11","CPSF2","LEO1","NUDT21","AHCYL1","CLP1","CSTF3","HSF1","SUPT5H","CTR9","NELFE","SCAF8","CPEB1","PAPOLB","SSU72P8","AP002990.1","SSU72P4","SSU72P5","SSU72P2","SSU72P7","SSU72P3","NCBP2","NCBP1","","CPSF4","NUDT21","CSTF2T","CPSF4")
# Transport of Mature Transcript to Cytoplasm	NUP160	TPR	THOC3	ZC3H11A	NDC1	U2AF2	NUP133	CPSF1	NUP37	THOC1	NUP50	AAAS	NUP188	POLDIP3	THOC5	SRSF5	RAE1	LUZP4	NUP93	CASC3	NUP88	MAGOHB	NUP107	SRSF9	SRSF3	NUP155	NCBP2	SRSF7	SRSF4	SRSF11	CPSF3	GLE1	NUP43	FYTTD1	DDX39A	SRSF6	NUP153	UPF3B	NUP85	THOC2	SYMPK	NUP214	THOC6	NUP210	NXT1	SRRM1	NUP42	SRSF1	WDR33	NCBP1	NUP54	DHX38	EIF4A3	FIP1L1	EIF4E	RANBP2	NUP205	SEC13	U2AF1	CHTOP	CPSF4	U2AF1L4	SRSF2	NXF1	MAGOH	NUP35	THOC7	SLBP	SLU7	CPSF2	CDC40	POM121	DDX39B	SARNP	RNPS1	NUP62	RBM8A	NXF2	NXF2B	POM121C
minor.spliceosome.go <- c("RNU4ATAC","RNU6ATAC","PHF5A","SF3B1","SF3B2","SF3B3","SF3B4","SF3B5","SF3B6","SNRPB","SNRPD1","SNRPD2","SNRPD3","SNRPE","SNRPF","SNRPG","SNRPN","DHX15","PDCD7","RNPC3","RNU11","RNU12","SNRNP25","SNRNP35","SNRNP48","YBX1","ZCRB1","ZMAT5","ZRSR2","CD2BP2","DDX23","EFTUD2","PRPF6","PRPF8","RNU5A-1","SNRNP200","SNRNP40","TXNL4A")
minor.splicing.go <- unique(c(gprofiler_database$`mRNA Splicing - Minor Pathway`, minor.spliceosome.go))
RNA.capping <-c("POLR2J","MNAT1","POLR2B","NCBP3","POLR2E","POLR2F","RNMT","POLR2C","ERCC2","POLR2I","GTF2H1","GTF2H3","RNGTT","NCBP2","GTF2F1","CDK7","CCNH","NCBP1","CMTR1","TGS1","POLR2D","GTF2H2","POLR2K","ERCC3","POLR2H","POLR2G","RAMAC","POLR2L","CMTR2","POLR2A","GTF2F2","SUPT5H","GTF2H4","RAMACL","GTF2H5")
exon.junction.complex <- c("UPF1","THRAP3","SMG6","TDRD3","PNN","R3HCC1","CASC3","MAGOHB","UPF3B","SRSF1","EIF4A3","SAP18","UPF2","MAGOH","R3HCC1L","UPF3A","PYM1","RNPS1","RBM8A","SRRM1","PRPF8","EIF4A3","MLN51","MAGOH","Y14","CASC3","EIF4A3","MAGOH","RBM8A","PYM1","NXF1","DDX39B")
RNA.polyadenylation <- c("PAF1","ZC3H3","PABPC1","CPSF1","SNRPA","YTHDC1","PAPOLA","ZC3H14","PABPC1L","CSTF2","RNF40","CPEB3","MTPAP","CPSF6","TENT4A","PAPOLG","CPSF3","TENT4B","SYMPK","CCNT1","GRSF1","CDC73","WDR33","CDK9","PNPT1","APP","FIP1L1","TUT1","CPSF7","RNF20","SSU72","TENT2","VIRMA","PCF11","CPSF2","LEO1","NUDT21","AHCYL1","CLP1","CSTF3","HSF1","SUPT5H","CTR9","NELFE","SCAF8","CPEB1","PAPOLB","SSU72P8","AP002990.1","SSU72P4","SSU72P5","SSU72P2","SSU72P7","SSU72P3","NCBP2","NCBP1","","CPSF4","NUDT21","CSTF2T","CPSF4")
exosome.complex <- c("MTREX","EXOSC7","EXOSC5","DIS3","KHSRP","GTPBP1","EXOSC3","EXOSC8","EXOSC9","EXOSC2","WDR74","ZFC3H1","MPHOSPH6","PNPT1","NVL","DIS3L2","CARHSP1","SUPV3L1","DIS3L","EXOSC1","EXOSC10","EXOSC4","C1D","EXOSC6")

phast <- phastCons100way.UCSC.hg38
gtf <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf.coords.RDS")

cleanIRFinder <- function(IRFinder = ac_who_vcp_vs_ctrl.IRFinder, deseq_res = NULL, mass_spec = NULL, species = "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,1], ensID = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,2], intron_type = str_split_fixed(IRFinder$Intron.GeneName.GeneID, "/", 3)[,3],
                        coords = paste(.$Chr, ":", .$Start, "-", .$End, sep = ""), Intron.GeneName.GeneID.Coords = paste(Intron.GeneName.GeneID, "/", coords, ":", Direction, sep = ""),
                        intron_length = GenomicRanges::width(as.GenomicRange(Intron.GeneName.GeneID.Coords)),
                        # reliable / stable intron expression
                        A.IR.coverage = A.SplicesMax + A.IntronDepth,
                        B.IR.coverage = B.SplicesMax + B.IntronDepth,
                        IR.coverage = (A.IR.coverage + B.IR.coverage)/2,
                        # reliable = case_when(A.IRok != "-" ~ "unreliable", B.IRok != "-" ~ "unreliable",  TRUE ~ "reliable"),
                        reliable = case_when(A.IRok %in% c("LowCover","LowSplicing") ~ "unreliable", B.IRok %in% c("LowCover","LowSplicing") ~ "unreliable",  TRUE ~ "reliable"),
                        # reliable = case_when(((A.SplicesMax + A.IntronDepth) > 100) & ((B.SplicesMax + B.IntronDepth) > 100) ~ "reliable", TRUE ~ "unreliable"), # based on https://github.com/lbroseus/SIRFindeR/blob/master/vignettes/SIRFindeR.pdf
                        # reliable_lenient = case_when(((A.SplicesMax + A.IntronDepth) > 10) & ((B.SplicesMax + B.IntronDepth) > 10) ~ "reliable", TRUE ~ "unreliable"),
                        retained = case_when((A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "retained", TRUE ~ "spliced"),
                        reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
                        # differential IR
                        ## effect size
                        IRratio.diff = A.IRratio - B.IRratio, # VCP - CTRL as per Dadi at IRFinder https://github.com/williamritchie/IRFinder/issues/99
                        IR.lfc = log(IRFinder$A.IRratio/(1-IRFinder$A.IRratio)) -  log(IRFinder$B.IRratio/(1-IRFinder$B.IRratio)), # as per lucile broseus https://github.com/lbroseus/SIRFindeR/issues/1#issuecomment-668611143
                        IRdirection = case_when(IRratio.diff > 0 ~ "up", TRUE ~ "down"), # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
                        ## significance
                        p0.05 = case_when(p.diff < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by p.diff < 0.05
                        # p0.05_IRr0.1 = case_when(p.diff < 0.05 & (A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and IRratio 0.1 in either sample
                        # p0.05_dIRr0.1 = case_when(p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
                        p0.05_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 ~ "IR up",
                                                      IRdirection == "down" & p.diff < 0.05 ~ "IR down",
                                                      TRUE ~ "none_significant"),
                        # p0.05_dIRr0.1_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR up",
                                                              # IRdirection == "down" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR down",
                                                              # TRUE ~ "none_significant"),
                        p0.05_reliable = case_when(p.diff < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
                        # p0.05_reliable_lenient = case_when(p.diff < 0.05 & reliable_lenient == "reliable" ~ "significant", TRUE ~ "none_significant"),
                        p0.05_reliable_IRdirection = case_when(IRdirection == "up" & p0.05_reliable == "significant" ~ "IR up",
                                                               IRdirection == "down" & p0.05_reliable  == "significant" ~ "IR down",
                                                               TRUE ~ "none_significant")) #,
                        # p0.05_reliable_lenient_IRdirection = case_when(IRdirection == "up" & p0.05_reliable_lenient == "significant" ~ "IR up",
                        #                                        IRdirection == "down" & p0.05_reliable_lenient  == "significant" ~ "IR down",
                        #                                        TRUE ~ "none_significant"),
                        # p0.05_reliable_IRr0.1 = case_when(p0.05_IRr0.1 == "significant" & p0.05_reliable == "significant" ~ "significant", TRUE ~ "none_significant")) # %>% left_join(gene_coordinates, by = "ensID") # takes too long
  gr_introns <- as.GenomicRange(IRFinder$Intron.GeneName.GeneID.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "Hsapiens"){
    print("GC content based on Hsapiens")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "Mmusculus"){
    print("GC content based on Mmusculus")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "Rnorvegicus"){
    print("GC content based on Rnorvegicus")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  if(species != "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = toupper(gene_name))  # covert gene_names to upper case for matching to Deseq results & Human gene names
  }
  # IRFinder$phast_cons <- gscores(phast, gr_introns) # too slow to keep in function
  if(!is.null(deseq_res)){
    IRFinder <- IRFinder %>% left_join(deseq_res, by = c("ensID", "gene_name")) # left join deseq2 results
    IRFinder <- IRFinder %>% dplyr::mutate(DGE.direction = case_when(log2FoldChange > 0 ~ "up", TRUE ~ "down"),
                                                           direction = case_when(log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
                                                                                 log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
                                                                                 log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
                                                                                 log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  }
  if(!is.null(mass_spec)){
    IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
    IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", vcp_vs_ctrl_ratio < 0 ~ "down"),
                                           massspec_ir_direction = case_when(vcp_vs_ctrl_ratio > 0 & IR.lfc < 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IR.lfc > 0 ~ "IR up & protein up",
                                                                          vcp_vs_ctrl_ratio < 0 & IR.lfc < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IR.lfc > 0 ~ "IR up & protein down"),
                                           massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
                                                                                        direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
                                                                                        direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
                                                                                        direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
                                                                                        TRUE ~ "other"))
  }
  return(IRFinder)
}


# clean the pooled replicate quantified IRFinder output (not differential IR). Needed to calculate absolute IR per condition.
cleanIRFinderQuant <- function(IRFinder = ac_who_ctrl.IRFinder, species = "Hsapiens"){
  IRFinder <- IRFinder %>% mutate(gene_name = str_split_fixed(IRFinder$Name, "/", 3)[,1], ensID = str_split_fixed(IRFinder$Name, "/", 3)[,2], intron_type = str_split_fixed(IRFinder$Name, "/", 3)[,3],
                                  coords = paste(.$Chr, ":", .$Start, "-", .$End, sep = ""), Name.Coords = paste(Name, "/", coords, ":", Strand, sep = ""),
                                  intron_length = GenomicRanges::width(as.GenomicRange(Name.Coords)),
                                  retained = case_when(IRratio >= 0.1 ~ "retained", TRUE ~ "spliced"),
                                  ExonDepth = case_when(SpliceRight > SpliceLeft ~ SpliceRight, TRUE ~ SpliceLeft), # aka SplicesMax
                                  # IR.coverage = ExonDepth + IntronDepth,
                                  # reliable = case_when(Warnings == "-" ~ "reliable", TRUE ~ "unreliable"), # case_when((ExonDepth + IntronDepth) > 100 ~ "reliable", TRUE ~ "unreliable"),
                                  reliable = case_when(Warnings %in% c("LowCover","LowSplicing") ~ "unreliable", TRUE ~ "reliable"), # case_when((ExonDepth + IntronDepth) > 100 ~ "reliable", TRUE ~ "unreliable"),
                                  # reliable_lenient = case_when((ExonDepth + IntronDepth) > 10 ~ "reliable", TRUE ~ "unreliable"),
                                  reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "reliable-retained", TRUE ~ "unreliable-or-spliced"))
                                  # reliable_lenient_retained = case_when(reliable_lenient == "reliable" & retained == "retained" ~ "reliable-retained", TRUE ~ "unreliable-or-spliced"))
  gr_introns <- as.GenomicRange(IRFinder$Name.Coords)
  seqlevelsStyle(gr_introns) <- "UCSC"
  if(species == "Hsapiens"){
    print("GC content based on Hsapiens")
    IRFinder$gc_content <- GCcontent(Hsapiens, gr_introns)
  }else if(species == "Mmusculus"){
    print("GC content based on Mmusculus")
    IRFinder$gc_content <- GCcontent(Mmusculus, gr_introns)
  }else if(species == "Rnorvegicus"){
    print("GC content based on Rnorvegicus")
    IRFinder$gc_content <- GCcontent(Rnorvegicus, gr_introns)
  }
  return(IRFinder)
}

# Compare the mutant and control pooled replicate quantified IRFinder output (not differential IR)
cleanIRFinderQuantCompared <- function(mutant = ac_who_vcp.IR, control = ac_who_ctrl.IR, differentialIR= NULL, deseq_res = NULL, mass_spec = NULL){
  IRFinder.control <- control %>% dplyr::select(ensID, gene_name, coords, Chr, Start, End, Name, Strand, intron_length, gc_content,
                                                "B.IRratio" = IRratio, "B.retained" = retained, "B.reliable" = reliable, "B.reliable_retained" = reliable_retained, "B.Coverage" = Coverage, "B.IntronDepth" = IntronDepth, 
                                                "B.SpliceLeft" = SpliceLeft, "B.SpliceRight" = SpliceRight, "B.SpliceExact" = SpliceExact, "B.Warnings" = Warnings, "B.intron_type" = intron_type)
  IRFinder.mutant <- mutant %>% dplyr::select(ensID, gene_name, coords, Chr, Start, End, Name, Strand, intron_length, gc_content,
                                              "A.IRratio" = IRratio, "A.retained" = retained, "A.reliable" = reliable, "A.reliable_retained" = reliable_retained, "A.Coverage" = Coverage, "A.IntronDepth" = IntronDepth, 
                                              "A.SpliceLeft" = SpliceLeft, "A.SpliceRight" = SpliceRight, "A.SpliceExact" = SpliceExact, "A.Warnings" = Warnings, "A.intron_type" = intron_type)
  IRFinder <- IRFinder.control %>% full_join(IRFinder.mutant, by = c("ensID", "gene_name", "coords", "Chr", "Start", "End", "Name", "Strand", "intron_length", "gc_content"))
  IRFinder <- IRFinder %>% mutate(
    # reliable event
    A.IR.coverage = max(A.SpliceLeft, A.SpliceRight) + A.IntronDepth,
    B.IR.coverage = max(B.SpliceLeft, B.SpliceRight) + B.IntronDepth,
    IR.coverage = (A.IR.coverage + B.IR.coverage)/2,
    reliable = case_when(A.Warnings %in% c("LowCover","LowSplicing") ~ "unreliable", B.Warnings %in% c("LowCover","LowSplicing") ~ "unreliable",  TRUE ~ "reliable"),
    retained = case_when((A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "retained", TRUE ~ "spliced"),
    reliable_retained = case_when(reliable == "reliable" & retained == "retained" ~ "retained", TRUE ~ "spliced"),
    # compare IR between conditions
    IRratio.diff = A.IRratio - B.IRratio, # Mutant - CTRL as per Dadi at IRFinder https://github.com/williamritchie/IRFinder/issues/99
    IR.lfc = log(IRFinder$A.IRratio/(1-IRFinder$A.IRratio)) -  log(IRFinder$B.IRratio/(1-IRFinder$B.IRratio)), # as per lucile broseus https://github.com/lbroseus/SIRFindeR/issues/1#issuecomment-668611143
    IR.direction = case_when(IRratio.diff > 0 ~ "up", TRUE ~ "down")) # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
  if(!is.null(differentialIR)){
    IRFinder.differential <- differentialIR %>% dplyr::select(Chr, Start, End, "Name" = Intron.GeneName.GeneID, "Strand" = Direction, p.diff)
    IRFinder <- IRFinder %>% left_join(IRFinder.differential, by = c("Chr", "Start", "End", "Name", "Strand")) # left join Differential IR results
    IRFinder <- IRFinder %>% dplyr::mutate(p0.05 = case_when(p.diff < 0.05 ~ "significant", TRUE ~ "none_significant"), # classify significant IR events by p.diff < 0.05
    # p0.05_IRr0.1 = case_when(p.diff < 0.05 & (A.IRratio >= 0.1 | B.IRratio >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and IRratio 0.1 in either sample
    # p0.05_dIRr0.1 = case_when(p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
    p0.05_IRdirection = case_when(IR.direction == "up" & p.diff < 0.05 ~ "IR up",
                                  IR.direction == "down" & p.diff < 0.05 ~ "IR down",
                                  TRUE ~ "none_significant"),
    # p0.05_dIRr0.1_IRdirection = case_when(IRdirection == "up" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR up",
    # IRdirection == "down" & p.diff < 0.05 & (abs(IRratio.diff) >= 0.1) ~ "IR down",
    # TRUE ~ "none_significant"),
    p0.05_reliable = case_when(p.diff < 0.05 & reliable == "reliable" ~ "significant", TRUE ~ "none_significant"),
    # p0.05_reliable_lenient = case_when(p.diff < 0.05 & reliable_lenient == "reliable" ~ "significant", TRUE ~ "none_significant"),
    p0.05_reliable_IR.direction = case_when(IR.direction == "up" & p0.05_reliable == "significant" ~ "IR up",
                                           IR.direction == "down" & p0.05_reliable  == "significant" ~ "IR down",
                                           TRUE ~ "none_significant"))
  }
  if(!is.null(deseq_res)){
    IRFinder <- IRFinder %>% left_join(deseq_res, by = c("ensID", "gene_name")) # left join deseq2 results
    IRFinder <- IRFinder %>% dplyr::mutate(DGE.direction = case_when(log2FoldChange > 0 ~ "up", TRUE ~ "down"),
                                           direction = case_when(log2FoldChange > 0 & IRratio.diff < 0 ~ "IR down & expression up",
                                                                 log2FoldChange > 0 & IRratio.diff > 0 ~ "IR up & expression up",
                                                                 log2FoldChange < 0 & IRratio.diff < 0 ~ "IR down & expression down",
                                                                 log2FoldChange < 0 & IRratio.diff > 0 ~ "IR up & expression down"))
  }
  if(!is.null(mass_spec)){
    IRFinder <- IRFinder %>%  left_join(mass_spec, by = c("gene_name"="name"))
    IRFinder <- IRFinder %>% dplyr::mutate(massspec_direction = case_when(vcp_vs_ctrl_ratio > 0 ~ "up", TRUE ~ "down"),
                                           massspec_ir_direction = case_when(vcp_vs_ctrl_ratio >= 0 & IR.lfc <= 0 ~ "IR down & protein up", vcp_vs_ctrl_ratio > 0 & IR.lfc > 0 ~ "IR up & protein up",
                                                                             vcp_vs_ctrl_ratio < 0 & IR.lfc < 0 ~ "IR down & protein down", vcp_vs_ctrl_ratio < 0 & IR.lfc > 0 ~ "IR up & protein down"),
                                           massspec_expression_ir_direction = case_when(direction == "IR down & expression up" & massspec_direction == "up" ~ "IR down & expression up & protein up",
                                                                                        direction == "IR down & expression up" & massspec_direction == "down" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "up" ~ "IR down & expression up & protein down", 
                                                                                        direction == "IR down & expression down" & massspec_direction == "down" ~ "IR down & expression down & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "up" ~ "IR up & expression down & protein up",
                                                                                        direction == "IR up & expression up" & massspec_direction == "up" ~ "IR up & expression up & protein up", 
                                                                                        direction == "IR up & expression up" & massspec_direction == "down" ~ "IR up & expression up & protein down",
                                                                                        direction == "IR up & expression down" & massspec_direction == "down" ~ "IR up & expression down & protein up",
                                                                                        TRUE ~ "other"))
  }
  return(IRFinder)
}


curated_profiles <- 
  function(gprofiles = IRevents_gprofiles_curated, colours = c("firebrick2", "dodgerblue3")){
    gprofiles %>% 
      #filter(term_name %in% curated_terms) %>%
      #separate(col = query, into = c("group", "direction"),sep ="_", remove = FALSE) %>%
      #filter( group == group_name )  %>% 
      distinct(term_name, .keep_all = TRUE) %>%
      arrange( -log10(p_value) ) %>% mutate(term_name = factor(term_name, levels = term_name)) %>%
      ggplot(aes( x = term_name, y = -log10(p_value), fill = query) ) + 
      geom_col(aes(y = ( (-log10(p_value))/sum(-log10(p_value)) * 100.0) )) + # makes more space to axis labels
      scale_fill_manual(values = colours) +
      # scale_fill_manual(values = get_palette("npg", 10)) +
      facet_grid(query ~ ., space = "free", scales = "free") +
      coord_flip() + theme_classic() + guides(fill=FALSE) +
      ylab(label = expression(log[10]~(P~value)) ) + xlab("") +
      theme(strip.background = element_blank(), legend.position = c(0.8, 0.2), strip.text.y = element_text(size=8, face = "bold",angle = 0), axis.text.y= element_text(size = 10)) +
      scale_y_continuous(expand = c(0, 0))
  }

curated_profiles_colours <- 
  function(gprofiles = IRevents_gprofiles_curated){
    gprofiles %>% 
      #filter(term_name %in% curated_terms) %>%
      #separate(col = query, into = c("group", "direction"),sep ="_", remove = FALSE) %>%
      #filter( group == group_name )  %>% 
      distinct(term_name, .keep_all = TRUE) %>%  arrange( -log10(p_value) ) %>%
      mutate(term_name = factor(term_name, levels = term_name)) %>%
      ggplot(aes( x = term_name, y = -log10(p_value), fill = query) ) + 
      geom_col(aes(y = ( (-log10(p_value))/sum(-log10(p_value)) * 100.0) )) + # makes more space to axis labels
      #scale_fill_manual(values = c("firebrick2", "dodgerblue")) +
      scale_fill_manual(values = get_palette("npg", 10)) +
      facet_grid(query ~ ., space = "free", scales = "free") +
      coord_flip() + theme_classic() + guides(fill=FALSE) +
      ylab(label = expression(log[10]~(P~value)) ) + xlab("") +
      theme(strip.background = element_blank()) +
      theme(strip.text.y = element_text(size=8, face = "bold",angle = 0)) +
      #labs(title = unique(group_name) ) +
      scale_y_continuous(expand = c(0, 0))
  }

labels = c("VIM", "COL1A1", "MYL6", "DDX5", "FN1", "FLNA", "HNRNPK", "FUS", "HNRNPL", "SFPQ", "EMP1", "HNRNPH1", "OSMR", "NONO", "HLA-A", "HLA-B", "HLA-C", "STAT3", "ANXA2", "SPARC", "HNRNPU", "HNRNPA2B1", "B2M", "PPIA", "LRP1","SRRM2", "SRSF5", "TNC", "SLC44A2", "CXCL2")



histogram_plot_irfinder <- 
  function(ctrl = ac_who_ctrl.IR, mutant = ac_who_vcp.IR, control = "ctrl", condition = "VCP"){
    ctrl.filt <- ctrl %>% dplyr::select(Chr, Start, End, Name, Strand, IRratio, gene_name, reliable) %>% mutate(condition = control)
    mutant.filt <- mutant %>% dplyr::select(Chr, Start, End, Name, Strand, IRratio, gene_name, reliable) %>% mutate(condition = condition)
    ctrl_and_mutant.filt <- ctrl.filt %>% bind_rows(mutant.filt)
    ggplot(filter(ctrl_and_mutant.filt, IRratio > 0), aes(x = IRratio, fill = reliable)) + 
      geom_histogram(colour="black", position="dodge",binwidth=0.03) +  theme_bw() + theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "right") + labs(fill = expression(Coverage)) + 
      geom_vline(xintercept = 0.1, linetype = 3, colour = "darkgrey")  + facet_grid(. ~ condition) + scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    } 



IRratio.diff_histogram_plot_irfinder <- 
  function(data = ac_who_vcp_vs_ctrl, reliabilityThr = 100){
    if(reliabilityThr == 100){
      ggplot(data, aes(x = IRratio.diff, fill = reliable)) + 
        geom_histogram(colour="black", position="dodge",binwidth=0.03) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top") + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 2, colour = "darkred") + geom_vline(xintercept = -0.1, linetype = 2, colour = "darkred")  + geom_vline(xintercept = 0, linetype = 1, colour = "black") +
        scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    } else{
      ggplot(data, aes(x = IRratio.diff, fill = reliable_lenient)) + 
        geom_histogram(colour="black", position="dodge",binwidth=0.01) +  theme_bw() + 
        theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank(), legend.position = "top") + labs(fill = expression(Coverage)) + 
        geom_vline(xintercept = 0.1, linetype = 2, colour = "darkred") + geom_vline(xintercept = -0.1, linetype = 2, colour = "darkred")  + geom_vline(xintercept = 0, linetype = 1, colour = "black")  +
        scale_fill_manual( values = c("dodgerblue", "firebrick2") )
    }
  }


reliabile_plot_irfinder <- 
  function(data = ac_who_ctrl, reliabilityThr = 100){
    ggplot(filter(data, IRratio > 0), aes( x = IRratio, y = (ExonDepth + IntronDepth),  label = gene_name)) + 
    theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank()) + theme_bw() +
    geom_point(colour = "purple", alpha = 0.5, size = 0.5) +  theme_bw() +  xlab( expression(IR~ratio)  ) +  ylab( expression(log[2]~Intron~Exon~abundance) ) + scale_y_continuous(trans='log2') +
    geom_hline(yintercept = log2(reliabilityThr), linetype = 3, colour = "darkred") + geom_vline(xintercept = 0.1, linetype = 3, colour = "darkred")
  }


ma_plot_irfinder <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = astrocyte.reactivity.labels){
    ggmaplot <- data %>% mutate(baseMean = log10(IR.coverage))
    labels_maplot <- data %>% mutate(baseMean = log10(IR.coverage)) %>% filter(p0.05_reliable == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = ggmaplot, aes(x = baseMean, y = IRratio.diff, colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        xlab(expression( Log[10]~coverage )) + 
        geom_text_repel(data = labels_maplot, aes(x = baseMean, y = IRratio.diff, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=3) +
        geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 1, linetype = 2)
  }

volcano_plot_irfinder <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = astrocyte.reactivity.labels, reliability = "reliable"){
    if(reliability == "reliable"){
      labels_plot <- data %>% filter(p0.05_reliable == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    } else if(reliability == "reliable_lenient"){
      labels_plot <- data %>% filter(p0.05_reliable_lenient == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_lenient_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
        scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
        ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
        geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0) +
        geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
    } else{
      labels_plot <- data %>% filter(p0.05 == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
      ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_reliable_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
      geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0) +
      geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }
  }

volcano_plot_irfinder_p0.05 <- 
  function(data = ac_ctrl_cyt_vs_nuc, xlabs = "Log[10] exon + intron coverage", labels_list = labels){
    labels_plot <- data %>% filter(p0.05 == "significant") %>% filter(gene_name %in% labels_list) %>% arrange( -(abs(IRratio.diff)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data = data, aes(x = IRratio.diff, y =  -log10(p.diff), colour = p0.05_IRdirection)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("firebrick2", "dodgerblue", "darkgray") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      ylab(expression( -log[10]~FDR )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
      geom_text_repel(data = labels_plot, aes(x = IRratio.diff, y = -log10(p.diff), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0) +
      geom_hline(yintercept =1.30, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }


scatter_plot_irfinder <-
  function(data = ctrl_vcp_scatter, labels_list = labels, x = IR.lfc_ctrl, y = IR.lfc_vcp){
    # labels_plot <- data %>% filter(gene_name %in% labels) %>% arrange( -(abs(x)) ) %>% distinct(gene_name, .keep_all = TRUE)
    ggplot(data, aes( x = x, y = y)) + geom_point(colour = "purple", alpha = 0.5, size = 0.5) +  theme_bw() +  
    # geom_text_repel(data = labels_plot, aes(x = x, y = y, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0.5) +
    theme(strip.background = element_rect(colour="black", fill="white"), panel.grid = element_blank()) +
    geom_hline(yintercept = 0, linetype = 3, colour = "darkgrey") + geom_vline(xintercept = 0, linetype = 3, colour = "darkgrey") + 
    geom_smooth(data = data, aes(x = x, y = y), method = "lm", se = FALSE, show.legend = TRUE, colour = "black") # + geom_abline(slope = -0.11, intercept = 0, linetype = 1)
  }



sinaplot <- 
  function(data = ac_who_ir_dge, groups = "p0.05_reliable_IRdirection", continuous = "log2FoldChange", stat.y.position = 1, tip.length = 0.01, stats.test = "t_test", x.position = NULL, dashed.line = 0, width = 0.2, flip = "no"){
    data2 <- data %>% dplyr::select(x = groups, y = continuous)
    print(paste("Stats test used = ", stats.test))
    if(stats.test == "t_test"){
      stats_test <- data2 %>% t_test(y ~ x, var.equal = TRUE) %>% adjust_pvalue() %>%  add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    } else{
      stats_test <- data2 %>% wilcox_test(y ~ x) %>% adjust_pvalue() %>%  add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    }
    if(flip == "yes"){
      if(dashed.line == "NA"){
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",4)) + scale_fill_manual(values = get_palette("npg",4)) +  
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8)) + coord_flip() +
          # geom_hline(yintercept = dashed.line, linetype = 2) + 
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position)
      } else{
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",4)) + scale_fill_manual(values = get_palette("npg",4)) +  
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8)) + coord_flip() +
          geom_hline(yintercept = dashed.line, linetype = 2) + coord_flip() +
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position) # add manually p values from binomial test
      }
    } else{
      if(dashed.line == "NA"){
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          # scale_colour_manual(values = get_palette("npg",4)) + scale_fill_manual(values = get_palette("npg",4)) +  
          scale_colour_manual(values = c("dodgerblue3", "firebrick2")) + scale_fill_manual(values = c("dodgerblue3", "firebrick2")) +  
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8)) + #coord_flip() +
          # geom_hline(yintercept = dashed.line, linetype = 2) + 
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position)
      } else{
        ggplot(data2, aes(x = x, y= y, color = x)) + 
          geom_violin() + geom_sina(size = 0.5) + geom_boxplot(data = data2, aes(fill = x), colour = "black", width=width, outlier.shape = NA) +
          # scale_colour_manual(values = get_palette("npg",4)) + scale_fill_manual(values = get_palette("npg",4)) +  
          scale_colour_manual(values = c("dodgerblue3", "firebrick2")) + scale_fill_manual(values = c("dodgerblue3", "firebrick2")) +  
          theme_classic() + theme(legend.position = "none", axis.text.x=element_text(size=10), axis.text.y=element_text(size=8)) + #coord_flip() +
          geom_hline(yintercept = dashed.line, linetype = 2) + #coord_flip() +
          stat_pvalue_manual(stats_test, label = "p.adj.signif", tip.length = tip.length, x = x.position) # add manually p values from binomial test
      }
    }    
  }

sinaplot.one.sample <- 
  function(data = boxplot.lfc.merged, groups = "category", continuous = "log2FoldChange", flip = "yes", tip.length = 0.01, labels = "gene_name", colours = 3, stat.y.position = 1.00, label_number = 10, stats.test = "t_test", expected = 0, width = 0.2, dashed.line = 0){
    data2 <- data %>% dplyr::select(group = groups, lfc = continuous, label = labels)
    print(data2)
    print(paste("Stats test used = ", stats.test))
    if(stats.test == "t_test"){
      stats_test <- data2 %>% dplyr::select(group, lfc) %>% group_by(group) %>% filter(lfc != Inf) %>% filter(lfc != -Inf) %>% t_test(lfc ~ 1, mu = expected) %>% add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    } else{
      stats_test <- data2 %>% dplyr::select(group, lfc) %>% group_by(group) %>% wilcox_test(lfc ~ 1, mu = expected) %>% add_significance() %>% mutate(y.position = stat.y.position)
      print(stats_test)
    }
    if(flip == "yes"){
      print(paste("Flipping coords"))
      if(colours < 4){
      ggplot(data2, aes(x = group, y= lfc)) + 
        geom_violin(aes(colour = group)) + 
        geom_sina(size = 0.5, aes(colour = group)) + 
        geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
        # scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
        scale_colour_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) + scale_fill_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) +
        theme_classic() + theme(legend.position = "none", axis.text.y=element_text(size=10)) +
        geom_hline(yintercept =dashed.line, linetype = 2) + coord_flip() +
        geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=12) +
        stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE)#, position_dodge(0.8)) # add manually p values
      } else{
        ggplot(data2, aes(x = group, y= lfc)) + 
          geom_violin(aes(colour = group)) + 
          geom_sina(size = 0.5, aes(colour = group)) + 
          geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
          # scale_colour_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) + scale_fill_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) +
          theme_classic() + theme(legend.position = "none", axis.text.y=element_text(size=10)) +
          geom_hline(yintercept =dashed.line, linetype = 2) + coord_flip() +
          geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=12) +
          stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE)#, position_dodge(0.8)) # add manually p values
      }
    } else{
      print(paste("Not flipping coords"))
      if(colours < 4){
      ggplot(data2, aes(x = group, y= lfc)) + 
        geom_violin(aes(colour = group)) + 
        geom_sina(size = 0.5, aes(colour = group)) + 
        geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
        # scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
        scale_colour_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) + scale_fill_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) +
        theme_classic() + theme(legend.position = "none", axis.text.y=element_text(size=10)) +
        geom_hline(yintercept = dashed.line, linetype = 2) + #coord_flip() +
        geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=12) +
        stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE)#, position_dodge(0.8)) # add manually p values
      } else{
        ggplot(data2, aes(x = group, y= lfc)) + 
          geom_violin(aes(colour = group)) + 
          geom_sina(size = 0.5, aes(colour = group)) + 
          geom_boxplot(data = data2, aes(fill = group), colour = "black", width=width, outlier.shape = NA) +
          scale_colour_manual(values = get_palette("npg",colours)) + scale_fill_manual(values = get_palette("npg",colours)) +
          # scale_colour_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) + scale_fill_manual(values = c("dodgerblue3", "firebrick2", "forestgreen")) +
          theme_classic() + theme(legend.position = "none", axis.text.y=element_text(size=10)) +
          geom_hline(yintercept = dashed.line, linetype = 2) + #coord_flip() +
          geom_text_repel(data = top_n(data2, label_number, abs(lfc)), aes(x = group, y = lfc, label = label), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=12) +
          stat_pvalue_manual(stats_test, label = "p.signif", tip.length = tip.length, xmin = "group", xmax = NULL, hide.ns = TRUE)#, position_dodge(0.8)) # add manually p values
      }
    }
  }
      
# Coverage Plots ---------------
options(ucscChromosomeNames=FALSE) # crucial to ensure AlignmentsTrack plot works
gtrack <- GenomeAxisTrack() # genome axis track
txdb <- makeTxDbFromGFF("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf", format="gtf") # make txdb from GTF - turn on only when needed
seqlevels(txdb, pruning.mode="coarse") <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y") # removes all the other chromsome annotations & scaffolds which are also in the BAM
ac_who_ctrl <- "/camp/home/ziffo/home/projects/astrocyte-meta-analysis/alignment/STAR/merged/ac_who_ctrl.bam"
ac_who_vcp <- "/camp/home/ziffo/home/projects/astrocyte-meta-analysis/alignment/STAR/merged/ac_who_vcp.bam"

getGeneIDsFromTxDb <- function(gr, txdb){
  stopifnot(is(gr, "GRanges"))
  stopifnot(length(gr)>0)
  stopifnot(is(txdb, "TxDb"))
  if(length(gr)>1){
    warning("The length of gr is greater than 1. Only first genomic location will be used.")
    gr <- gr[1]
  }
  genes <- genes(txdb, columns="gene_id")
  genes <- subsetByOverlaps(genes, gr)
  return(genes$gene_id)
}

# Track viewer https://www.bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#introduction
track_viewer <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, xpad = 1000){
  top_event <- txt %>% dplyr::filter(p0.05 == "significant") %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  IRratio_ctrl <- top_event %>% pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  # intron coords
  thechr <-  top_event %>% pull(Chr) # without chr before number
  st <-  top_event %>% pull(Start) # 35184593
  en <-  top_event %>% pull(End) # 35187000
  strand <- top_event %>% pull(Direction) 
  # whole gene coords
  # gene_test_gtf <- gtf %>% filter(gene_name %in% gene_test)
  # gene_chr <- gene_test_gtf %>% pull(seqnames) %>% unique
  # gene_st <- gene_test_gtf %>% pull(start) %>% min
  # gene_en <- gene_test_gtf %>% pull(end) %>% max
  # gr <- GRanges(gene_chr, IRanges(gene_st, gene_en), strand= strand) # set width of tracks as entire gene
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  # trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% pull(ensID) # get ENSG of gene_test
  gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  # View the tracks
  optSty <- optimizeStyle(trackList(VCP, CTRL, gene_track))
  trackList <- optSty$tracks
  names(trackList)[3] <- gene_test # rename gene_track with gene_test name - shows as the gene label
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .1, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .4)
  setTrackStyleParam(trackList[[2]], "height", .4)
  setTrackStyleParam(trackList[[3]], "height", .1)
  setTrackStyleParam(trackList[[3]], "ylabgp", list(cex=1.5))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[3]], "color", c("dodgerblue3", "black"))
  
  # save figure
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder/figures/coverage/", gene_test, "_coverage_plot.png", sep = ""), height = 3.5, width = 10, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.545, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.545, y=.2, just="bottom")
  dev.off()
}

track_viewer_tx <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, xpad = 1000, path = NULL){
  txt %>% dplyr::filter(p0.05 == "significant" & gene_name == gene_test) %>% print
  top_event <- txt %>% dplyr::filter(p0.05 == "significant" & gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff))
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable IR events. Only the top deltaIR event will be plotted.", sep = ""))
  }
  IRratio_ctrl <- top_event %>% pull(B.IRratio) # B = ctrl
  IRratio_ctrl <- format(round(IRratio_ctrl,2),nsmall = 2)
  IRratio_vcp <- top_event %>% pull(A.IRratio) # A = vcp
  IRratio_vcp <- format(round(IRratio_vcp,2),nsmall = 2)
  thechr <-  top_event %>% pull(Chr) # without chr before number
  st <-  top_event %>% pull(Start) # 35184593
  en <-  top_event %>% pull(End) # 35187000
  strand <- top_event %>% pull(Direction) 
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand) # Granges for entire gene (not just the retained intron)
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% pull(ensID) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  if (length(trackList) > 8 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue3", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.6, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_coverage_plot.png", sep = ""), height = 2.2, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 2, lty = "dashed")
  grid.text(IRratio_ctrl, x=.57, y=.65, just="bottom")
  grid.text(IRratio_vcp, x=.57, y=.2, just="bottom")
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}

# Plot multiple IR events on a single gene coverage plot
track_viewer_tx_multiple <- function(gene_test = "SFPQ", txt = ac_who_vcp_vs_ctrl, mut.bam = ac_who_vcp, ctrl.bam = ac_who_ctrl, path = NULL){
  txt_nrow <- txt  %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable IR events.", sep = ""))
  }
  events <- txt %>% dplyr::filter(p0.05_reliable == "significant" & gene_name == gene_test) 
  
  # set GRanges for entire gene
  # whole gene coords
  gene_test_gtf <- gtf %>% filter(gene_name %in% gene_test)
  gene_chr <- gene_test_gtf %>% pull(seqnames) %>% unique %>% as.numeric
  gene_st <- gene_test_gtf %>% pull(start) %>% min %>% as.numeric
  gene_en <- gene_test_gtf %>% pull(end) %>% max %>% as.numeric
  strand <- gene_test_gtf %>% pull(strand) %>% unique %>% as.character
  gr <- GRanges(gene_chr, IRanges(gene_st, gene_en), strand= strand) # set width of tracks as entire gene
  CTRL <- importBam(ctrl.bam, ranges=gr, pairs = TRUE) # import BAM
  VCP <- importBam(mut.bam, ranges=gr, pairs = TRUE) # import BAM
  CTRL$dat <- coverageGR(CTRL$dat) # calculate coverage
  VCP$dat <- coverageGR(VCP$dat) # calculate coverage
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(VCP, CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  # setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  # setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .1, 0.03, .01))
  # setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  # setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .15, -0.1, .01))
  setTrackStyleParam(trackList[[1]], "height", .36)
  setTrackStyleParam(trackList[[2]], "height", .36)
  setTrackXscaleParam(trackList[[1]], "draw", TRUE)
  setTrackXscaleParam(trackList[[1]], "gp", list(cex=0.9))
  setTrackStyleParam(trackList[[1]], "color", c("darkgrey", "black"))
  setTrackStyleParam(trackList[[2]], "color", c("darkgrey", "black"))
  
  # only show first 6 transcript tracks - remove remaining transcript tracks
  if (length(trackList) > 8 ){
    print("removing transcript tracks > 6:")
    for(i in 9:length(trackList)) {
      print(trackList[[i]]$name)
    }
    idx = c(9:length(trackList))
    trackList <- trackList[-idx]
    # [[9:length(trackList)]] <- NULL
  }
  # modify the transcript tracks
  for(i in 3:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue3", "black"))
  }
  
  print("Extracting each IR event coords.")
  event <- list()
  IRratio_ctrl <- list()
  IRratio_vcp <- list()
  thechr <- list()
  st <- list()
  en <- list()
  strand <- list()
  for (i in 1:nrow(events)) {
    event[[i]] <- events[i, ] # iterate through rows of events df  
    IRratio_ctrl[[i]] <- event[[i]] %>% pull(B.IRratio) # B = ctrl
    IRratio_ctrl[[i]] <- format(round(IRratio_ctrl[[i]],2),nsmall = 2)
    IRratio_vcp[[i]] <- event[[i]] %>% pull(A.IRratio) # A = vcp
    IRratio_vcp[[i]] <- format(round(IRratio_vcp[[i]],2),nsmall = 2)
    thechr[[i]] <-  event[[i]] %>% pull(Chr) # without chr before number
    st[[i]] <-  event[[i]] %>% pull(Start) # 35184593
    en[[i]] <-  event[[i]] %>% pull(End) # 35187000
    strand[[i]] <- event[[i]] %>% pull(Direction) 
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  for (i in 1:nrow(events)) {
    addGuideLine(c(st[[i]], en[[i]]), vp=vp, col = "firebrick2", lwd = 1, lty = "dashed")
    # grid.text(IRratio_ctrl[[i]], x=.545, y=.6, just="bottom")
    # grid.text(IRratio_vcp[[i]], x=.545, y=.2, just="bottom")
  }
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  
  print("Saving figure.")
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder/figures/coverage/", path, gene_test, "_tx_multiple_coverage_plot.png", sep = ""), height = 2.5, width = 8, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  for (i in 1:nrow(events)) {
    addGuideLine(c(st[[i]], en[[i]]), vp=vp, col = "firebrick2", lwd = 1, lty = "dashed")
    # grid.text(IRratio_ctrl[[i]], x=.545, y=.6, just="bottom")
    # grid.text(IRratio_vcp[[i]], x=.545, y=.2, just="bottom")
  }
  grid.text(gene_test, x=.09, y=.87, just="bottom", gp=gpar(cex=1.5))
  dev.off()
}


track_viewer_fractions <- function(gene_test = "RND3", txt = ac_vcp_vs_ctrl.filt, mut.nuc.bam = ac_nuc_vcp, ctrl.nuc.bam = ac_nuc_ctrl, mut.cyt.bam = ac_cyt_vcp, ctrl.cyt.bam = ac_cyt_ctrl, xpad = 1000, path = NULL){
  txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% print
  top_event <- txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% top_n(., 1, abs(IRratio.diff.nuc))
  txt_nrow <- txt %>% drop_na %>% arrange(IRratio.diff.nuc) %>% dplyr::filter(gene_name == gene_test) %>% nrow
  if(txt_nrow > 1){
    print(paste(gene_test, " has ", txt_nrow, " significant reliable nuclear IR events. Only the top deltaIR event will be plotted.", sep = ""))
  }
  IRratio_nuc_ctrl <- top_event %>% pull(IRratio.ctrl.nuc) # B = ctrl
  IRratio_nuc_ctrl <- format(round(IRratio_nuc_ctrl,2),nsmall = 2)
  IRratio_nuc_vcp <- top_event %>% pull(IRratio.vcp.nuc) # A = vcp
  IRratio_nuc_vcp <- format(round(IRratio_nuc_vcp,2),nsmall = 2)
  IRratio_cyt_ctrl <- top_event %>% pull(IRratio.ctrl.cyt) # B = ctrl
  IRratio_cyt_ctrl <- format(round(IRratio_cyt_ctrl,2),nsmall = 2)
  IRratio_cyt_vcp <- top_event %>% pull(IRratio.vcp.cyt) # A = vcp
  IRratio_cyt_vcp <- format(round(IRratio_cyt_vcp,2),nsmall = 2)
  
  thechr <-  top_event %>% pull(Chr) # without chr before number
  st <-  top_event %>% pull(Start) # 35184593
  en <-  top_event %>% pull(End) # 35187000
  strand <- top_event %>% pull(Direction) 
  
  gr <- GRanges(thechr, IRanges(st - xpad, en + xpad), strand= strand) # Granges for entire gene (not just the retained intron)
  # import BAM
  nuc.CTRL <- importBam(ctrl.nuc.bam, ranges=gr, pairs = TRUE)
  nuc.VCP <- importBam(mut.nuc.bam, ranges=gr, pairs = TRUE)
  cyt.CTRL <- importBam(ctrl.cyt.bam, ranges=gr, pairs = TRUE)
  cyt.VCP <- importBam(mut.cyt.bam, ranges=gr, pairs = TRUE)
  # calculate coverage
  nuc.CTRL$dat <- coverageGR(nuc.CTRL$dat) 
  nuc.VCP$dat <- coverageGR(nuc.VCP$dat)
  cyt.CTRL$dat <- coverageGR(cyt.CTRL$dat)
  cyt.VCP$dat <- coverageGR(cyt.VCP$dat)
  # Build Gene model
  trs <- geneModelFromTxdb(txdb, org.Hs.eg.db, gr = gr) # build gene model showing transcripts
  # ids <- getGeneIDsFromTxDb(gr, txdb) # get ENSG of gene_test
  ids <- top_event %>% pull(ensID) # get ENSG of gene_test
  # gene_track <- geneTrack(ids,txdb)[[1]] # gene all exons on single line
  
  print("Creating the coverage and gene model tracks.")
  optSty <- optimizeStyle(trackList(cyt.VCP, cyt.CTRL, nuc.VCP, nuc.CTRL, trs))
  trackList <- optSty$tracks
  viewerStyle <- optSty$style
  setTrackViewerStyleParam(viewerStyle, "xaxis", FALSE)
  setTrackViewerStyleParam(viewerStyle, "margin", c(.01, .10, -0.1, .01))
  # optimise coverage plots
  for(i in 1:4) {
    # setTrackStyleParam(trackList[[i]], "draw", TRUE)
    setTrackXscaleParam(trackList[[i]], "gp", list(cex=0.9))
    setTrackStyleParam(trackList[[i]], "height", .18)
    setTrackStyleParam(trackList[[i]], "color", c("darkgrey", "black"))
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0.8))
  }
  # optimise the transcript tracks
  if (length(trackList) > 10 ){
    print(paste("There are ", length(trackList), " ensembl transcripts. Plotting only the first 6.")) 
    for(i in 11:length(trackList)) {
      print(trackList[[i]]$name) 
    }
    idx = c(11:length(trackList))
    trackList <- trackList[-idx] # only plot first 6 transcript tracks - remove remaining transcript tracks
  }
  for(i in 5:length(trackList)) {
    setTrackStyleParam(trackList[[i]], "height", .03)
    setTrackStyleParam(trackList[[i]], "ylabgp", list(cex=0))
    setTrackStyleParam(trackList[[i]], "color", c("dodgerblue3", "black"))
  }
  
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_cyt_vcp, x=.545, y=.15, just="bottom")
  grid.text(IRratio_cyt_ctrl, x=.545, y=.35, just="bottom")
  grid.text(IRratio_nuc_vcp, x=.545, y=.55, just="bottom")
  grid.text(IRratio_nuc_ctrl, x=.545, y=.75, just="bottom")
  grid.text(gene_test, x=.05, y=.87, just="bottom", gp=gpar(cex=1.0))
  
  print("Saving figure.")
  png(file = paste("/camp/home/ziffo/home/projects/astrocyte-meta-analysis/splicing/IRFinder/figures/coverage/fractions/", path, gene_test, "_tx_coverage_plot.png", sep = ""), height = 4, width = 7, units ="in", res = 300)
  vp <- viewTracks(trackList, gr=gr, viewerStyle=viewerStyle)
  addGuideLine(c(st, en), vp=vp, col = "firebrick2", lwd = 3, lty = "dashed")
  grid.text(IRratio_cyt_vcp, x=.545, y=.15, just="bottom")
  grid.text(IRratio_cyt_ctrl, x=.545, y=.35, just="bottom")
  grid.text(IRratio_nuc_vcp, x=.545, y=.55, just="bottom")
  grid.text(IRratio_nuc_ctrl, x=.545, y=.75, just="bottom")
  grid.text(gene_test, x=.06, y=.87, just="bottom", gp=gpar(cex=1.0))
  dev.off()
}

# HTSeq annotation
# mus_musculus.gtf <- import("/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gtf") %>% as.tibble
# mus.musculus.gene2ens <- mus_musculus.gtf %>% dplyr::select(gene_id, gene_name)
mus.musculus.gene2ens <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Mus_musculus.GRCm38.101.gene2ens.RDS")
# homo_sapiens.gtf <- import("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gtf") %>% as.tibble
# gene2ens <- homo_sapiens.gtf %>% dplyr::select(gene_id, gene_name)
gene2ens <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.gene2ens.RDS") %>% unique
trans2gene <- readRDS("/camp/home/ziffo/home/genomes/ensembl/Homo_sapiens.GRCh38.99.trans2ens.RDS") %>% unique

ens_to_gene <- function(ens) {
  gene <- replace(ens, ens %in% as.character(gene2ens$gene_id), as.character(gene2ens$gene_name[gene2ens$gene_id %in% ens]))
  # gene <- gene2ens %>% dplyr::filter(gene_id == ens) %>% pull(gene_name) %>% as.character()
  return(gene)
}
gene_to_ens <- function(gene) {
  ens <- gene2ens %>% dplyr::filter(gene_name == gene) %>% pull(gene_id) %>% as.character()
  print(ens)
}
ens_to_gene_list <- function(ens) {
  gene <- list()
  for (i in seq_along(ens)) {
    gene[[i]] <- replace(ens[[i]], ens[[i]] %in% as.character(gene2ens$gene_id), as.character(gene2ens$gene_name[gene2ens$gene_id %in% ens[[i]]]))
  }
  return(gene)
}


cleanVASTTOOLS <- function(diff = ac_who_vcp_vs_ctrl_diff, comp = ac_who_vcp_vs_ctrl_compare){
  txt <- diff %>% left_join(dplyr::select(comp, -starts_with("nuclear_"),-starts_with("cytoplasmic_")), by = c("GENE", "EVENT"))
  txt <- txt %>% mutate(type = case_when(grepl("HsaEX", EVENT) ~ "exon skipping",
                                         grepl("HsaINT", EVENT) ~ "intron retention",
                                         grepl("HsaALTA", EVENT) ~ "alternate 3'",
                                         grepl("HsaALTD", EVENT) ~ "alternate 5'"))
  txt <- txt %>% mutate(compare.direction = case_when(dPSI > 0 ~ "up", TRUE ~ "down"),
                        dPSI0.1 = case_when(dPSI > 0.1 ~ "significant", TRUE ~ "none_significant"), 
                        dPSI0.1_direction = case_when(compare.direction == "up" & dPSI0.1 > 0.1 ~ "up",
                                                      compare.direction == "down" & dPSI0.1 < -0.1 ~ "down",
                                                      TRUE ~ "none_significant"),
                        diff.direction = case_when(E.dPsi. > 0 ~ "up", TRUE ~ "down"), # if positive then IRratio greater in VCP, negative then IRratio greater in CTRL
                        MV0.1 = case_when(MV.dPsi._at_0.95 > 0.1 ~ "significant", TRUE ~ "none_significant"), 
                        MV0.1_dPSI0.1 = case_when(MV.dPsi._at_0.95 > 0.1 & abs(E.dPsi.) >= 0.1 ~ "significant", TRUE ~ "none_significant"), # significant IR events p.diff < 0.05 and delta IRratio of 0.1 between samples
                        MV0.1_direction = case_when(diff.direction == "up" & MV.dPsi._at_0.95 > 0.1 ~ "up",
                                                    diff.direction == "down" & MV.dPsi._at_0.95 > 0.1 ~ "down",
                                                    TRUE ~ "none_significant"),
                        MV0.1_dPSI0.1_direction = case_when(diff.direction == "up" & MV.dPsi._at_0.95 > 0.1 & (abs(E.dPsi.) >= 0.1) ~ "up",
                                                            diff.direction == "down" & MV.dPsi._at_0.95 > 0.1 & (abs(E.dPsi.) >= 0.1) ~ "down",
                                                            TRUE ~ "none_significant"))
  fromSplit <- str_split_fixed(txt$COORD, ":", 2)
  txt$Chr <- fromSplit[,1] %>% gsub("chr", "", .)
  coord <- str_split_fixed(fromSplit[,2], "-", 2)
  txt$Start <- coord[,1] %>% as.numeric
  txt$End <- coord[,2] %>% as.numeric
  return(txt)
}


makePieChart <- function(txt = ac_who_vcp_vs_ctrl.VASTTOOLS){
  num.events <- nrow(txt)
  print(paste0("total number of events: ", num.events ))
  # treat each class of event as separate
  txt <- table(unlist(txt$type) ) 
  print(txt)
  print( sum(txt) )
  txt.plot <- c() 
  txt.plot["exon skipping"] <- sum( c(txt["exon skipping"]), na.rm = TRUE )
  txt.plot["alternate 3'"] <- sum( c(txt["alternate 3'"]), na.rm = TRUE )
  txt.plot["alternate 5'"] <- sum( c( txt["alternate 5'"] ), na.rm = TRUE )
  txt.plot["intron retention"] <- sum( c(txt["intron retention"]), na.rm = TRUE)
  txt.plot <- as.data.frame(txt.plot) 
  names(txt.plot) <- "varcounts" 
  print( txt.plot)
  txt.plot$type <- paste(row.names(txt.plot), " n = ", txt.plot$varcounts, sep = "")
  txt.plot <- txt.plot[ order(txt.plot$varcounts,decreasing=FALSE),]  
  txt.plot$type <- factor(txt.plot$type, levels = rev(txt.plot$type) )
  txt.plot <- dplyr::mutate(txt.plot, pos = cumsum(varcounts) - 0.5*varcounts) 
  txt.plot$prop <- signif( (txt.plot$varcounts / sum(txt.plot$varcounts) ) * 100, 3) 
  txt.plot$prop <- format(round(txt.plot$prop, 1), nsmall = 1)
  txt.plot$prop <- paste0(txt.plot$prop, "%")
  txt.plot <- txt.plot[ txt.plot$varcounts >0 ,]
  txt.plot$label <- paste(txt.plot$type, "\n (", txt.plot$prop, ")", sep = "")
  txt.plot$label <- factor(txt.plot$label, levels = rev(txt.plot$label) )
  txt.plot$label <- gsub("n = ", "\n n = ", txt.plot$label)
  pie <- ggplot(txt.plot, aes(x="", y = varcounts, fill = type)) + 
    geom_bar(width = 1, stat="identity")  + 
    geom_text(aes(x=1.9, y = pos, label = label), size = 4 ) +
    coord_polar(theta="y") +
    #scale_fill_brewer("",palette="Paired", direction = 1) + 
    #scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Paired"))(13)) + 
    scale_fill_manual(values = get_palette("npg", 11)) + scale_color_manual(values = get_palette("npg", 11)) +
    theme_void() + theme(legend.position = "none") # + ggtitle(title) + ggplot2::annotate("text", x = 1.8, y = sum(txt.plot$varcounts) / 2, label = paste0("total number of events at padj < ", FDRlimit,": ", num.events ) )
  return(pie)
}


nest_batch <- function(inner, within, set_to=".") {
  n_instances <- apply(table(inner, within)!=0, 1, sum)
  if (any(n_instances>1)) {
    stop(paste(names(n_instances)[n_instances>1], collapse=", "), " appear in multiple parents")
  }
  levels_in <- levels(inner)
  # each inner level has in a unique group, find it
  corresponding_group <-within[match(levels_in, inner)]
  # for each group, find index of first inner level
  which_inner <- match(levels(within), corresponding_group)
  # and set it to the common value
  if (any(levels(inner)[-which_inner]==set_to )) { # we're going to duplicate an existing level so best stop
    stop("One of the batches is already called ", set_to, ". Please use a non-existing level")
  }
  levels(inner)[which_inner] <- set_to
  inner
}

res_top <- function(r, i=1) {
  ind <- which(r$padj<0.05) # indexes of significant genes
  ind[order(abs(r$log2FoldChange)[ind])][i]
}


ma_plot_deseq <- 
  function(data = mn_d25.res, labels_list = ma.labels){
    labels_maplot <- data %>% filter(gene_name %in% ma.labels)
    ggplot(data = data, aes(x = baseMean, y = log2FoldChange, colour = padj < 0.05)) +  geom_point(size = 0.5) +
      scale_colour_manual( values = c("darkgray", "firebrick2") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      xlab(expression( Log[2]~mean~expression )) + ylab(expression(Log[2]~fold~change~VCP:CTRL)) + 
      geom_text_repel(data = labels_maplot, aes(x = baseMean, y = log2FoldChange, label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=3) +
      geom_hline(yintercept = 0, linetype = 2)
  }


volcano_plot_deseq <- 
  function(data = mn_d25.res, labels_list = volcano.labels){
    labels_plot <- data %>% filter(gene_name %in% volcano.labels)
    ggplot(data = data, aes(x = log2FoldChange, y =  -log10(pvalue), colour = padj < 0.05)) +  geom_point(size = 0.5) + #xlim(1,5) + ylim(-3,3) +
      scale_colour_manual( values = c("darkgray", "firebrick2") ) +  theme_bw() +  theme(panel.grid = element_blank(), legend.title = element_blank(), legend.position = "none") +
      ylab(expression( -log[10]~P~value )) + #xlab( expression( IR~log[2]~fold~change~CTRL~cytoplasmic:nuclear)) +
      geom_text_repel(data = labels_plot, aes(x = log2FoldChange, y = -log10(pvalue), label = gene_name), fontface = "italic", colour = "black", min.segment.length = 0, box.padding = 0.5, point.padding = 0, size=3) +
      geom_hline(yintercept =5, linetype = 2) + geom_vline(xintercept =0, linetype = 2)
  }

make_matrix <- function(df,rownames = NULL){
  my_matrix <-  as.matrix(df)
  if(!is.null(rownames))
    rownames(my_matrix) = rownames
  my_matrix
}


left_join_NA <- function(x, y, ...) {
  left_join(x = x, y = y, by = ...) %>% 
    mutate_each(funs(replace(., which(is.na(.)), 0)))
}

