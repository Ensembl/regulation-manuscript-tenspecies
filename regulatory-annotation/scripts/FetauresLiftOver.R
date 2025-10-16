# Script: FeaturesLiftOver.R
# Author: Gabriela Merino - Ensembl Regulation
# Date: 10-2025
# Description: Project ensembl regulatory features between two assemblies.

suppressPackageStartupMessages({
  library(dplyr)
  library(argparse)
  library(stringr)
  library(bedtoolsr)
  library(tidyr)
  set.seed(1234)
})
options(scipen = 999)

parsingArgs <- function(){
  parser <- ArgumentParser(prog = 'Rscript FeaturesLiftOver.R',
                           description ='Lifting over regulatory features 
                           between two assemblies')
  parser$add_argument('fromAnnotation', metavar = 'fromAnnotation', 
                      type = "character", default = NULL, nargs = 1,
                      help = 'Full path to the gff3 file containing the 
                      regulatory features to be lifted over.')
  parser$add_argument('fromMotifs', metavar = 'fromMotifs', 
                      type = "character", default = NULL, nargs = 1,
                      help = 'Full path to the gff3 file containing the TF 
                      motifs features to be lifted over.')
  parser$add_argument('toAnnotation', metavar = 'toAnnotation', 
                      type = "character", default = NULL, nargs = 1,
                      help = 'Name of the gff3 file where lifted 
                      regulatory features will be exported.')
  parser$add_argument('toMotifs', metavar = 'toMotifs', 
                      type = "character", default = NULL, nargs = 1,
                      help = 'Name of the gff3 file where lifted 
                      TF motifs features will be exported.')
  parser$add_argument('chainFile', metavar = 'chainFile', 
                      type = "character", nargs = 1, 
                      help = 'Full path to the chain file.')
  parser$add_argument('outPath', metavar = 'outPath', 
                      type = "character", nargs = 1, 
                      help = 'Full path to the output folder.')
  parser$add_argument('fromAssemblyInfo', metavar = 'fromAssemblyInfo', 
                      type = "character", nargs = 1,  
                      help = 'Full path to the report file of the original
                      assembly used for RFs identification, useful for
                      chromosomes ID mapping')
  parser$add_argument('toAssemblyInfo', metavar = 'toAssemblyInfo', 
                      type = "character", nargs = 1,  
                      help = 'Full path to the report file of the new
                      assembly, useful for chromosomes ID mapping')
  parser$add_argument('annotatedTSSs', metavar = 'TSSsFile', type='character',
                      nargs = 1, help = 'Full path to the bed file 
                      containing TSSs (chr, start, end, strand) and gene 
                      annotations of the new assembly.')
  parser$add_argument('toAssemblyChr', metavar = 'toAssemblyChr', 
                      type = "character", nargs = 1,  
                      help = 'Full path to the chromsomes file of the new
                      assembly.')
  parser$add_argument('--keepRFsInDiffChr', metavar = 'keepRFsInDiffChr', 
                      action = 'store', default = FALSE,
                      help = 'Should features mapping to different 
                      chromosomes be kept? Default is to remove them.')
  parser$add_argument('--updateGenePromLink', metavar = 'updateGenePromLink', 
                      action = 'store', default = TRUE,
                      help = 'Should the gene-promoter links be updated? 
                      Default is yes, keep the links that are still valid in
                      the new assembly. If it is FALSE, then new links will
                      be created'.)  
  parser$add_argument('--keepNonTSSsProms', metavar = 'keepNonTSSsProms', 
                      action = 'store', default = FALSE,
                      help = 'Should promoters do not overlapping TSSs 
                      (without genes linked to) be kept? Default is FALSE')
  parser$add_argument('--promColor', metavar = 'promColor', 
                      action = 'store', default = '#ff0000',
                      help = 'Color of promoters.')
  return(parser)
}


loadRFs <- function(annotationFile, assemblyInfo = NULL, 
                    chrIDType = NULL, extractPromoters = FALSE){
  RFs <- read.delim(annotationFile, skip = 1, header = F, sep = '\t')
  if(extractPromoters){
    # promoters extended regions will also be lifted over    
    RFs %>%
      bind_rows(RFs %>% 
                  filter(V3 == 'promoter') %>%
                  mutate(V3 = 'extended_promoter',
                         V4 = as.numeric(str_split(str_split(V9,
                                                             'extended_start=',
                                                             simplify = T)[,2],
                                                   ';', simplify = T)[,1]),
                         V5 = as.numeric(str_split(str_split(V9, 
                                                             'extended_end=', 
                                                             simplify = T)[,2],
                                                   ';', simplify = T)[,1]))
      )-> RFs
  }
  if(!is.null(assemblyInfo) & !is.null(chrIDType)){
    refName <- names(assemblyInfo)[1]
    assemblyInfo <- assemblyInfo[assemblyInfo[, refName] %in% RFs[, 1], ]
    RFs[, 1] <- factor(RFs[, 1],
                       levels = assemblyInfo[, 1])
    levels(RFs[,1]) <- assemblyInfo[, chrIDType]
    
  }
  return(RFs)
}
removeRFsFromDifferentChr <- function(rfs, fromRFs){
  # remove features that were lifted over to a different chromosome
  fromRFs %>% 
    rename(fromChr = V1) %>%
    select(ID, fromChr) -> fromRFsInfo
  rfs %>%
    left_join(fromRFsInfo) %>%
    filter(V1 == fromChr) %>%
    select(-fromChr) -> rfs
}
removePromotersWithoutTSS <- function(rfs, annotatedTSSsAndInfo){
  rfs %>%
    filter(V3 == 'promoter') %>%
    select(V1, V4, V5, ID) %>%
    mutate(V4 = V4 - 1) %>% 
    # identify those promoters that do not overlap a TSS
    bt.intersect(annotatedTSSsAndInfo, wa = T, v = T) -> promsToRemove
  if(nrow(promsToRemove) >0){
    promsToRemove %>%
      pull(V4) -> promsToRemove
  }else{promsToRemove<- ''}
  rfs %>% 
    filter(!ID %in% promsToRemove) -> rfs
  return(rfs)
}
updatingGenePromLink <-function(rfs, annotatedTSSsAndInfo, 
                                promColor = promColor){
  # flagging ensembl gene ids linked to the same gene name
  annotatedTSSsAndInfo %>%
    group_by(V1, V7) %>%
    mutate(duplicated = if_else(length(unique(V5)) > 1, TRUE,
                                FALSE)) -> TSSAndInfo
  
  # extracting gene annotation from the original rfs file
  rfs %>%
    filter(V3 == 'promoter' & grepl('gene_name', V9)) %>%
    mutate(gene_id = str_split(str_split(V9, 'gene_id=', simplify = T)[,2],
                               ';', simplify = T)[,1],
           gene_name = str_split(str_split(V9, 'gene_name=', simplify = T)[,2],
                                 ';', simplify = T)[,1],
           gene_biotype = str_split(str_split(V9, 'gene_biotype=', simplify = T)[,2],
                                    ';', simplify = T)[,1]) %>%
    bind_rows(rfs %>% 
                filter(V3 == 'promoter' & !grepl('gene_name', V9)) %>%
                mutate(gene_id = str_split(str_split(V9, 'gene_id=',
                                                     simplify = T)[,2],
                                           ';', simplify = T)[,1],
                       gene_name = '',
                       gene_biotype = str_split(str_split(V9, 'gene_biotype=',
                                                          simplify = T)[,2],
                                                ';', simplify = T)[,1]))->proms
  proms %>%
    select(V1, gene_id, gene_name, gene_biotype) %>%
    rename(chr = V1) %>%
    separate_longer_delim(c(gene_id, gene_name, gene_biotype), delim=',') %>%
    unique() %>% 
    as.data.frame() -> geneInfoFrom
  
  # extract overlapping between promoters and TSSs in the new annotation
  proms %>% 
    mutate(V4 = V4 -1,
           extended_start = as.numeric(str_split(str_split(V9, 
                                                           'extended_start=',
                                                           simplify = T)[,2],
                                                 ';', simplify = T)[,1]),
           extended_end = as.numeric(str_split(str_split(V9, 
                                                         'extended_end=',
                                                         simplify = T)[,2],
                                               ';', simplify = T)[,1])) %>%
    select(V1, V4, V5, extended_start, extended_end, ID, gene_id,
           gene_name, gene_biotype) %>%
    separate_longer_delim(c(gene_id, gene_name, gene_biotype), delim=',') %>%
    bt.intersect(TSSAndInfo[,c(-4)], wa = TRUE, 
                 wb = TRUE) -> proms
  # remove links that are absent in the new assembly
  proms %>%  
    filter(((V7 == V13 | V8 == V15) & !V16) |
             (V7 == V13 & V8 == V15))%>%
    mutate(V15 = if_else(V13 == V15, '', V15)) %>%
    mutate(V7 = if_else(V7 == V13, V7, '')) %>%
    mutate(V8 = if_else(V8 == V15, V8, '')) %>%
    mutate(V9 = if_else(V9 == V14, V9, '')) %>%
    select(V6, V4, V5, V7, V8, V9) %>%
    unique() %>%
    rename(ID = V6, extended_start = V4, extended_end = V5, gene_id = V7,
           gene_name = V8, gene_biotype = V9) %>%
    group_by(ID, extended_start, extended_end) %>%
    summarise(gene_id = paste(gene_id, collapse =','),
              gene_name = paste(gene_name, collapse =','),
              gene_biotype = paste(gene_biotype, collapse =','), 
              .groups = 'drop') %>%
    mutate(attr = paste(paste('ID=', ID, sep =''), 
                        paste('extended_start=', extended_start, sep =''),
                        paste('extended_end=', extended_end, sep =''), 
                        sep = ';')) %>%
    mutate(attr = if_else(str_remove_all(gene_id,',') == '',
                          attr, paste(attr, ';gene_id=', gene_id, sep = ''))) %>%
    mutate(attr = if_else(str_remove_all(gene_name,',') == '',
                          attr, paste(attr, ';gene_name=', 
                                      gene_name, sep = ''))) %>%
    mutate(attr = if_else(str_remove_all(gene_biotype,',') == '',
                          attr, paste(attr, ';gene_biotype=', gene_biotype,
                                      sep = ''))) %>%
    mutate(attr = paste(attr, paste('color=', promColor, sep = ''),
                        sep = ';')) %>%
    select(ID, attr) %>%
    as.data.frame() -> proms
  proms %>%
    left_join(rfs %>% 
                filter(V3 == 'promoter')) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, attr, ID) %>%
    rename(V9 = attr) -> proms
  rfs %>%
    filter(V3 != 'promoter') %>%
    bind_rows(proms)-> rfs
  return(rfs)
}
addingGenePromLink <-function(rfs, annotatedTSSsAndInfo, promColor = promColor){
  # adding gene info for promoters
  rfs %>% 
    filter(V3 == 'promoter') %>%
    mutate(extended_start = as.numeric(str_split(str_split(V9, 
                                                           'extended_start=',
                                                           simplify = T)[,2],
                                                 ';', simplify = T)[,1]) - 1,
           extended_end = as.numeric(str_split(str_split(V9, 
                                                         'extended_end=',
                                                         simplify = T)[,2],
                                               ';', simplify = T)[,1])) %>%
    select(V1, V4, V5, extended_start, extended_end, ID) %>%
    bt.intersect(annotatedTSSsAndInfo, wa = TRUE, wb = TRUE) -> proms
  if(ncol(proms) == 13){
    proms %>% mutate(V13 = if_else(V13 == V11, '', V13)) -> proms
  }else{
    proms %>%
      bind_cols('V13' = '') -> proms
  }
  proms %>%  
    select(c(V4, V5, V6, V11, V12, V13)) %>%
    unique() %>%
    group_by(V6, V4, V5) %>%
    summarise(gene_id = paste(V11, collapse =','),
              gene_name = paste(V13, collapse =','),
              gene_biotype = paste(V12, collapse =','), .groups = 'drop') %>%
    mutate(attr = paste(paste('ID=', V6, sep =''), 
                        paste('extended_start=', V4 + 1 , sep =''),
                        paste('extended_end=', V5, sep =''), sep = ';')) %>%
    mutate(attr = paste(attr, 
                        if_else(str_remove_all(gene_name,',') == '',
                                paste(';gene_id=', gene_id, sep = ''),
                                paste(paste(';gene_id=', gene_id, sep = ''),
                                      paste('gene_name=', gene_name, sep = ''), 
                                      sep = ';')), sep ='')) %>%
    mutate(attr = paste(attr, paste('gene_biotype=', gene_biotype, sep = ''), 
                        paste('color=', promColor, sep = ''), sep = ';')) %>%
    select(V6, attr) %>%
    as.data.frame() -> proms
  rfs %>%
    left_join(proms %>% rename(ID = V6)) %>%
    mutate(V9 = if_else(!is.na(attr), attr, V9)) %>%
    select(-attr) -> rfs
  return(rfs)
}

removeOverlappingRFs <-function(rfs, fromRFs, tfms, fromTFMs){
  # promoters, open chromatin regions and enhancer should not overlap between 
  # them. 
  # if no new links have been added, the promoters links have been updated to
  # match the current annotation so we don't expect any overlap.
  # For EMARs, we discard the ones overlapping removed promoters. We will
  # pick the one overlapping at least one promoter in the original annotation,
  # otherwise the longest one
  # For enhancers and TF motifs, we remove the ones overlapping removed EMARs.
  # For CTCFs, we remove the ones that do not overlap any kept CTCF motif.
  rfs %>%
    filter(V3 == 'promoter') %>%
    mutate(extended_start = as.numeric(str_split(str_split(V9, 
                                                           'extended_start=',
                                                           simplify = T)[,2],
                                                 ';', simplify = T)[,1]),
           extended_end = as.numeric(str_split(str_split(V9, 
                                                         'extended_end=',
                                                         simplify = T)[,2],
                                               ';', simplify = T)[,1])) %>%
    rename(chr = V1) %>%
    select(chr, extended_start, extended_end, ID) -> promotersToKeep
  # Checking overlapping, we keep promoters that overlap protein_coding gene  
  # TSSs over lncRNA
  promotersToKeep %>%
    bt.sort() %>%
    bt.merge() -> promMerg
  if(nrow(promMerg) < nrow(promotersToKeep)){
  # overlapping to fix
    promMerg %>%
      bt.intersect((promotersToKeep %>%
                      left_join(rfs %>% select(ID, V9)) %>%
                      mutate(biotype = str_split(str_split(V9, 'gene_biotype=', 
                                                           simplify = T)[,2],
                                                 ';', simplify = T)[,1]) %>%
                      select(chr, extended_start, extended_end, ID, biotype)), 
                   wa = T, wb = T) %>%
      group_by( V1, V2, V3) %>%
      mutate(n =n()) %>% filter(n>1) %>% 
      mutate(keep = grepl('protein_coding', V8) | 
               !any(grepl('protein_coding', V8))) %>% 
      filter(!keep) %>%
      pull(V7) -> promsLncRNAToRemove 
    promotersToKeep %>%
      filter(!ID %in% promsLncRNAToRemove) -> promotersToKeep
  }
  # we identify the EMARs overlapping the promoters in the original assembly
  fromRFs %>% 
    filter(V3 == 'EMAR') %>%
    select(V1, V4, V5, ID) %>%
    bt.intersect((fromRFs %>% 
                    filter(ID %in% (promotersToKeep %>% pull(ID))) %>% 
                    select(V1, V4, V5, ID)), wa = T, u = T) %>%
    pull(V4) -> EMARsOvPromToKeep
  
  # We check EMARs overlapping, if it exists, we try to keep the one overlapping
  # a promoter linked to a protein_coding gene. Otherwise, we keep the one 
  # one overlapping the higher number of promoters, or the one having the longer
  # region overlapping promoters
  rfs %>% 
    filter(V3 == 'EMAR' & ID %in% EMARsOvPromToKeep) %>%
    select(V1, V4, V5) %>%
    bt.merge() %>%
    bt.intersect((rfs %>% 
                    filter(ID %in% EMARsOvPromToKeep) %>% 
                    select(V1, V4, V5, ID)), wa = T, wb = T) %>%
    group_by(V1, V2, V3) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    mutate(width=V6-V5) %>%
    select(V1,V2,V3,width, V7) %>%
    rename(chr=V1, start=V2,end=V3,ID=V7)-> EMARsOvToCheck
  # now we overlap them with the promoters in the old coordinates
  fromRFs %>% 
    filter(ID %in% (EMARsOvToCheck %>% pull(ID))) %>%
    select(V1, V4, V5, ID) %>%
    bt.intersect((fromRFs %>% 
                    filter(ID %in% (promotersToKeep %>% pull(ID))) %>%
                    mutate(biotype = str_split(str_split(V9, 'gene_biotype=', 
                                                         simplify = T)[,2],
                                               ';', simplify = T)[,1]) %>%
                    select(V1, V4, V5, ID, biotype)), wb = T) %>%
    select(V2, V3, V4, V8, V9) %>%
    mutate(width = V3-V2) %>%
    group_by(V4) %>%
    mutate(nProm = n(), 
              width = sum(width),
              nPC = sum(grepl('protein_coding', V9)))%>%
    rename(ID=V4, promID = V8) %>%
    select(ID, nProm, width, nPC, promID) %>%
    left_join(EMARsOvToCheck %>% select(chr, start, end, ID)) %>%
    arrange(chr, start, end, desc(nPC), desc(nProm), desc(width)) %>%
    group_by(chr,start,end) %>%
    mutate(keep = (ID == ID[1])) %>%
    filter(!keep) %>%
    rename(EMARID = ID) %>% 
    select(EMARID, promID) -> EMARsAndPromOvBetweenThemToRemove
  
  promotersToKeep %>%
    filter(!ID %in% (EMARsAndPromOvBetweenThemToRemove %>% 
                       pull(promID))) -> promotersToKeep
  
  # we need to remove EMARs overlapping kept promoters in the translated 
  # coordinates but that originally do not overlap them
  rfs %>% 
    filter(V3 == 'EMAR') %>%
    select(V1, V4, V5, ID) %>%
    bt.intersect((rfs %>% 
                    filter(ID %in% (promotersToKeep %>% pull(ID))) %>% 
                    select(V1, V4, V5, ID)), wa = T, u = T) %>%
    filter(!V4 %in% EMARsOvPromToKeep) %>%
    select(V4)  %>%
    bind_rows(EMARsAndPromOvBetweenThemToRemove %>%
                rename(V4=EMARID) %>% select(V4)) %>%
    pull(V4) %>%
    unique() -> EMARsOvPromToRemove
  # we keep those EMARs or the first one when no overlapping with a kept
  # promoter is present
  rfs %>%
    filter(V3 == 'EMAR' & !(ID %in% EMARsOvPromToRemove)) %>%
    select(V1, V4, V5, ID) %>%
    bt.merge() %>%
    bt.intersect((rfs %>% 
                    filter(V3 == 'EMAR'& !(ID %in% EMARsOvPromToRemove)) %>% 
                    select(V1, V4, V5, ID)), wa = T, wb = T) %>%
    mutate(ovPromoter = V7 %in% EMARsOvPromToKeep) %>%
    arrange(V1, V2, V3, V7) %>%
    group_by(V1, V2, V3) %>%
    summarise(V7 = ifelse(any(ovPromoter),
                          V7[ovPromoter][1],
                          V7[1]), .groups = 'drop') %>%
    as.data.frame() %>%
    rename(ID = V7) %>% select(ID) %>%
    left_join(rfs) %>% 
    select(V1, V4, V5, ID)-> EMARsToKeep
  
  fromRFs %>%
    filter(V3 == 'EMAR' & !(ID %in% 
                              (EMARsToKeep %>% pull(ID)))) %>%
    select(V1, V4, V5, ID) -> EMARsToRemove

  # Because some promoters may not overlap an EMAR in the original assembly or
  # the EMAR could not be lifted over, we need to check that we are removing 
  # those that are now overlapping a kept EMAR in the new one.
  promotersToKeep %>%
    bt.intersect(rfs %>%
                   filter(V3 == 'EMAR' &
                            !(ID %in% (EMARsToRemove %>% pull(ID)))) %>%
                   select(V1, V4, V5, ID), wa = T, wb = T) %>%
    group_by(V1, V2, V3, V4) %>%
    summarise(EMAR=paste(V8, collapse=';')) %>%
    rename(ID = V4) %>%
    select(V1,V2, V3, ID, EMAR) %>%
    left_join(fromRFs %>%
                filter(ID %in% (promotersToKeep %>% pull(ID))) %>%
                select(V1,V4,V5,ID) %>%
                bt.intersect(fromRFs %>%
                               filter(V3 == 'EMAR') %>%
                               select(V1, V4, V5, ID), wa = T, wb = T) %>%
                select(V4, V8) %>%
                rename(ID=V4, fromEMAR=V8)) %>%
    mutate(keep = str_detect(EMAR, fromEMAR)) %>%
    filter(!keep) %>%
    pull(ID) -> promotersToRemove
  
  promotersToKeep %>%
    filter(!ID %in% promotersToRemove ) -> promotersToKeep
  

  fromRFs %>%
    filter(V3 == 'enhancer') %>%
    select(V1, V4, V5, ID) %>%
    bt.intersect(EMARsToRemove, wa = T, u = T) %>%
    pull(V4) -> enhancersToRemove

  rfs %>%
    filter(V3 == 'enhancer' & !(ID %in% enhancersToRemove)) %>%
    select(V1, V4, V5, ID) -> enhancerToKeep#-> enh
  
  rfs %>%
    filter(V3 == 'open_chromatin_region') %>%
    select(V1, V4, V5, ID) %>%
    # remove OCRs overlapping promoters
    bt.intersect(promotersToKeep %>% select(chr, extended_start,
                                            extended_end),
                 wa = T, v = T) %>%
    # remove OCRs overlapping enhancers
    bt.intersect(enhancerToKeep, wa = T, v = T) -> ocrs
  ocrs %>%
    bt.merge() %>%
    bt.intersect(ocrs, wa = T, wb = T) %>%
    arrange(V1, V2, V3, V4, V5, V6) %>%
    group_by(V1, V2, V3) %>%
    summarise(V7 = V7[1], .groups = 'drop') %>%
    as.data.frame() %>%
    rename(ID = V7) %>% select(ID) -> OCRToKeep
  
  rfs %>%
    filter(V3 != 'CTCF_binding_site') %>%
    left_join(promotersToKeep %>%
                select(ID) %>%
                bind_rows((enhancerToKeep %>% select(ID))) %>%
                bind_rows((OCRToKeep %>% select(ID))) %>%
                bind_rows((EMARsToKeep %>% select(ID))) %>%
                mutate(status='keep')) %>%
    filter(!is.na(status)) %>%
    select(-status) -> keptRFs
  
  # we will remove CTCFs and TFMs that overlap regions whose features have been 
  # removed
  fromRFs %>%
    filter(V3 != 'CTCF_binding_site') %>%
    left_join(promotersToKeep %>%
                select(ID) %>%
                bind_rows((enhancerToKeep %>% select(ID))) %>%
                bind_rows((OCRToKeep %>% select(ID))) %>%
                bind_rows((EMARsToKeep %>% select(ID))) %>%
                mutate(status='keep')) %>%
    filter(is.na(status)) %>%
    select(V1, V4, V5, V3, ID) %>%
    bt.intersect((fromRFs %>%
                    filter(V3 != 'CTCF_binding_site') %>%
                    left_join(promotersToKeep %>%
                                select(ID) %>%
                                bind_rows((enhancerToKeep %>% select(ID))) %>%
                                bind_rows((OCRToKeep %>% select(ID))) %>%
                                bind_rows((EMARsToKeep %>% select(ID))) %>%
                                mutate(status='keep')) %>%
                    filter(!is.na(status)) %>%
                    select(V1, V4, V5)), wa = T, v = T) -> removedRFsFrom
  
  fromRFs %>% 
    filter(V3 == 'CTCF_binding_site') %>%
    select(V1, V4, V5, V7, ID) %>%
    bt.intersect(removedRFsFrom, wa = T, u = T) %>%
    pull(V5) -> ctcfsToRemove 
  
  fromTFMs %>% 
    select(V1, V4, V5, V7, ID) %>%
    bt.intersect(removedRFsFrom, wa = T, u = T) %>%
    pull(V5) -> TFMsToRemove  
  # CTCFs should not overlap between them. 
  rfs %>%
    filter(V3 == 'CTCF_binding_site' & !(ID %in% ctcfsToRemove)) %>%
    select(V1, V4, V5, V7, ID) %>%
    bt.intersect((rfs %>% 
                    filter(V3 == 'CTCF_binding_site' & !(ID %in% ctcfsToRemove)) %>% 
                    select(V1, V4, V5, V7, ID)), wa = T, wb = T) -> ctcfs2consider
  
  # TFMs predicted by using the same PWM should show a perfect overlap.
  # We need to randomly pick one of these motifs. 
  tfms %>%
    mutate(PWM = str_split(str_split(V9, 'matrix_id=', simplify = T)[,2],
                           ';', simplify = T)[,1]) -> tfms
  tfms %>%
    filter(!(ID %in% TFMsToRemove)) %>% 
    select(V1, V4, V5, ID, PWM, V7) %>%
    bt.intersect((tfms %>% 
                    filter(!(ID %in% TFMsToRemove)) %>% 
                    select(V1, V4, V5, ID, PWM, V7)), s = T, wa = T,
                 wb = T) -> tfms2consider
  tfms2consider %>%
    group_by(V1, V2, V3, V4, V5, V6, V8, V9, V11) %>%
    summarise(V7 = unique(V7),
              V10 = V10[1], .groups = 'drop') %>%
    filter(V1 == V7, V2 == V8, V3 == V9, V5 == V11) %>%
    pull(V10) %>%
    unique() -> tfms2consider
  # we remove CTCFs do not overlapping kept TFMs
  fromRFs %>%
    filter(V3 == 'CTCF_binding_site' & !(ID %in% ctcfsToRemove)) %>%
    select(V1, V4, V5, ID) %>%
    bt.intersect(fromTFMs %>%
                   filter(ID %in% tfms2consider &
                            (grepl('CTCF', V9) | grepl('Ctcf', V9))) %>%
                   select(V1, V4, V5), wa = T, u = T) %>%
    pull(V4) -> ctcfsIntfms
  rfs %>%
    filter(ID %in% ctcfsIntfms) %>%
    bind_rows(keptRFs) %>% 
    arrange(V1,V4,V5)-> finalLORFs
  tfms %>%
    filter(ID %in% tfms2consider) %>%
    arrange(V1,V4,V5, V7) %>% 
    select(-PWM) -> finalLOTFMs
  return(list(finalLORFs = finalLORFs, finalLOTFMs = finalLOTFMs))
}

main <- function(parser){
  args <- parser$parse_args()
  fromAnnotation <- args$fromAnnotation
  toAnnotation <- args$toAnnotation
  fromMotifs <- args$fromMotifs
  toMotifs <- args$toMotifs
  outPath <- args$outPath
  if(!file.exists(outPath)){
    system(paste('mkdir ', outPath, sep=''))
  }
  
  fromAssemblyInfo <- args$fromAssemblyInfo
  toAssemblyInfo <- args$toAssemblyInfo
  toAssemblyChr <- args$toAssemblyChr
  
  promColor <- args$promColor
  # information to update chromosomes names to match the chain file
  conn <- file(fromAssemblyInfo, open="r")
  linn <- readLines(conn)
  nl <- 1
  headerFound <- FALSE
  assemblyMapHeader <- NULL
  while(!headerFound){
    if(grepl('#', str_split(linn[nl], ' ', simplify = T)[,1])){
      if(grepl('\t', linn[nl])){
        assemblyMapHeader <- str_split(str_split(linn[nl], ' ', 
                                                 simplify = T)[,-1], 
                                       '\t', simplify = T)
      }
      nl <- nl + 1
    }else{
      headerFound <- TRUE
    }
  }
  close(conn)
  fromAssemblyMap <- read.table(fromAssemblyInfo, sep = '\t', 
                                comment.char = "#")
  colnames(fromAssemblyMap) <- assemblyMapHeader
  annotatedTSSs <- args$annotatedTSSs
  chainFile <- args$chainFile
  keepRFsInDiffChr <- as.logical(args$keepRFsInDiffChr)
  keepNonTSSsProms <- as.logical(args$keepNonTSSsProms)
  updateGenePromLink <- as.logical(args$updateGenePromLink)
  
  
  toChr <- tryCatch(read.delim(toAssemblyChr, header=F, sep= "\t")[,1], 
                                   warning = function(w) NULL) 
  stopifnot("File listing the chromosomes names of the new assembly doesn't exist" = 
              !is.null(toChr))
  
  annotatedTSSsAndInfo <- tryCatch(read.delim(annotatedTSSs, header=F, 
                                              sep= "\t"), 
                                   warning = function(w) NULL) 
  stopifnot("File with annotated TSSs is required to update promoter-gene links" = 
              !is.null(annotatedTSSsAndInfo))
  
  # Inferring chromosome IDs styles based on the chain file
  incon <- gzfile(chainFile)
  aux <- str_split(readLines(incon,1), ' ', simplify = TRUE)[3]
  chrIDTypeInLiftOver <- names(fromAssemblyMap[grepl(aux, fromAssemblyMap)])
  
  #load the annotation file
  fromRFs <- loadRFs(fromAnnotation)
  fromRFs %>% 
    mutate(ID = as.character(str_split(str_split(V9, 'ID=', 
                                                 simplify = T)[,2],
                                       ';', simplify = T)[,1])) ->fromRFs
  fromTFMs <- loadRFs(fromMotifs)
  fromTFMs %>% 
    mutate(ID = as.character(str_split(str_split(V9, 'ID=', 
                                                 simplify = T)[,2],
                                       ';', simplify = T)[,1])) -> fromTFMs
  # load the annotation file, add extended regions for promoters and update 
  # chromosomes names 
  oldRFs <- loadRFs(fromAnnotation, fromAssemblyMap, chrIDTypeInLiftOver, 
                      extractPromoters = TRUE)
  oldTFMs <- loadRFs(fromMotifs, fromAssemblyMap, chrIDTypeInLiftOver)

  write.table(oldRFs, paste(outPath, 'tmp.gff', sep = '/'), sep = '\t',
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  # liftOver the annotation
  
  system(paste('liftOver -gff', paste(outPath, 'tmp.gff', sep = '/'),
               chainFile, paste(outPath, toAnnotation, sep = '/'),
               paste(outPath, 'tmpLiftOver', sep ='/'), collapse = " "))
  write.table(oldTFMs, paste(outPath, 'tmp.gff', sep = '/'), sep = '\t',
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  system(paste('liftOver -gff', 
               paste(outPath, 'tmp.gff', sep = '/'), chainFile, 
               paste(outPath, toMotifs, sep = '/'), 
               paste(outPath, 'tmpLiftOver', sep ='/'), collapse = " "))
  
  conn <- file(toAssemblyInfo, open="r")
  linn <- readLines(conn)
  nl <- 1
  headerFound <- FALSE
  assemblyMapHeader <- NULL
  while(!headerFound){
    if(grepl('#', str_split(linn[nl], ' ', simplify = T)[,1])){
      if(grepl('\t', linn[nl])){
        assemblyMapHeader <- str_split(str_split(linn[nl], ' ', 
                                                 simplify = T)[,-1], 
                                       '\t', simplify = T)
      }
      nl <- nl + 1
    }else{
      headerFound <- TRUE
    }
  }
  close(conn)
  toAssemblyMap <- read.table(toAssemblyInfo, sep = '\t', comment.char = "#")
  colnames(toAssemblyMap) <- assemblyMapHeader
  aux <- str_split(readLines(incon,1), ' ', simplify = TRUE)[8]
  chrIDTypeOutLiftOver <- names(toAssemblyMap[grepl(aux, toAssemblyMap)])

  toAssemblyMap <- toAssemblyMap[, c(chrIDTypeOutLiftOver,
                                      names(toAssemblyMap)[1])]
  newRFs <- loadRFs(paste(outPath, toAnnotation, sep = '/'), 
                    assemblyInfo = toAssemblyMap, 
                    chrIDType = names(toAssemblyMap)[2],
                    extractPromoters = FALSE)
  newTFMs <- loadRFs(paste(outPath, toMotifs, sep = '/'), 
                    assemblyInfo = toAssemblyMap, 
                    chrIDType = names(toAssemblyMap)[2],
                    extractPromoters = FALSE)
  # keep features in canonical chromosomes
  newRFs %>%
    filter( V1 %in% toChr) %>%
    mutate(ID = str_split(str_split(V9, 'ID=', simplify = T)[,2],
                          ';', simplify = T)[,1]) -> newRFs
  newTFMs %>%
    filter( V1 %in% toChr) %>%
    mutate(ID = str_split(str_split(V9, 'ID=', simplify = T)[,2],
                          ';', simplify = T)[,1]) -> newTFMs
  # update promoters core and extended regions
  newRFs %>% 
    filter(grepl('promoter', V3)) %>% 
    group_by(V1, V2, V6, V7, V8, V9, ID) %>%
    summarise(extended_start = min(V4), extended_end = max(V5),
              V4 = max(V4), V5 = min(V5), .groups = 'drop') %>%
    as.data.frame() %>%
    mutate(V3 = 'promoter',
           oldExtEnd = str_split(str_split(V9, 'extended_end=', simplify=T)[,2],
                                 ';', simplify = T)[,1]) %>%
    mutate(gene_link_info = str_split(V9, paste('extended_end=', oldExtEnd,
                                                sep =''), 
                                      simplify = T)[,2]) %>%
    mutate(V9 = paste('ID=', ID, ';extended_start=', extended_start,
                      ';extended_end=', extended_end, gene_link_info,
                      sep = '')) %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, ID) %>%
    bind_rows(newRFs %>% filter(!grepl('promoter', V3))) %>%
    arrange(V1, V4, V5) -> loRFs
  
  # if specified, features mapping to different chromosomes will be removed
  if(!keepRFsInDiffChr){
    loRFs <- removeRFsFromDifferentChr(loRFs, fromRFs)
    newTFMs <- removeRFsFromDifferentChr(newTFMs, fromTFMs)
  }
  # if specified, promoters do not overlapping TSSs will be removed
  if(!keepNonTSSsProms){
    loRFs <- removePromotersWithoutTSS(loRFs, annotatedTSSsAndInfo)
  }
    
  # if specified, add new promoter-gene links
  if(updateGenePromLink){
    loRFs <- updatingGenePromLink(loRFs, annotatedTSSsAndInfo, promColor)
  }else{
    loRFs <- addingGenePromLink(loRFs, annotatedTSSsAndInfo, promColor)  
  }
  #remove overlapping features
  finalSetOfFeatures <- removeOverlappingRFs(loRFs, fromRFs, newTFMs, fromTFMs)
  finalLORFs <- finalSetOfFeatures[['finalLORFs']]
  finalLOTFMs <- finalSetOfFeatures[['finalLOTFMs']]
  
  finalLORFs[,'V1'] <- factor(finalLORFs[,'V1'], levels = toChr) 
  finalLOTFMs[,'V1'] <- factor(finalLOTFMs[,'V1'], levels = toChr)
  finalLORFs %>% 
    arrange(V1, V4, V5) %>%
    select(-ID) -> finalLORFs
  finalLOTFMs %>% 
    arrange(V1, V4, V5) %>%
    select(-ID) -> finalLOTFMs
  # save final gff3 file
  gff3h <- "##gff-version 3"
  RFFile <- paste(outPath, toAnnotation, sep = '/')
  system(paste('rm', RFFile, sep =' '))
  RFFileCon <- file(RFFile)
  writeLines(text = gff3h, con = RFFileCon)
  close(RFFileCon)
  write.table(finalLORFs, file = RFFile, 
              col.names = F, row.names = F, sep = "\t", quote = F,
              append = T)
  TFMFile <- paste(outPath, toMotifs, sep = '/')
  system(paste('rm', TFMFile, paste(outPath, 'tmp.gff', sep = '/'),
               paste(outPath, 'tmpLiftOver', sep ='/'), sep =' '))
  TFMFileCon <- file(TFMFile)
  writeLines(text = gff3h, con = TFMFileCon)
  close(TFMFileCon)
  write.table(finalLOTFMs, file = TFMFile, 
              col.names = F, row.names = F, sep = "\t", quote = F,
              append = T)
}
parser <- parsingArgs()
main(parser)