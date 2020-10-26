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

