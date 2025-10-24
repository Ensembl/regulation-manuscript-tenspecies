# Script: RegulatoryBuild.R
# Author: Gabriela Merino - Ensembl Regulation
# Date: 10-2025
# Description: Regulatory features identification with the Ensembl regulatory build.

suppressPackageStartupMessages({
  library(dplyr)
  library(bedtoolsr)
  library(argparse)
  library(stringr)
})
options(scipen = 999)
parsingArgs <- function(){
  parser <- ArgumentParser(prog = 'Rscript RegulatoryBuild.R',
                           description = 'Regulatory features annotation')
  
  parser$add_argument('species', metavar = 'S', type = "character", nargs =1,
                      help = 'Species name in lower case and using underscore
                      as a separator (e.g Xxx_yyy).')
  parser$add_argument('assembly', metavar = 'A', type = "character", nargs =1,
                      help = 'Assembly name.')
  parser$add_argument('ESpePre', metavar = 'EPrefix', type = "character", 
                      nargs = 1, help = 'Ensembl prefix for the species.')
  parser$add_argument('startID', metavar = 'startID', type = "character", 
                      nargs = 1,
                      help = 'Integer specifying the starting value for
                      feature IDs.')
  parser$add_argument('annotatedTSSs', metavar = 'TSSsFile', type = "character", 
                      nargs = 1, help = 'Full path to the bed file containing 
                      TSSs annotation (chr, start, end, strand).')
  parser$add_argument('annotatedCDS', metavar = 'CDSFile', type = "character", 
                      nargs = 1, help = 'Full path to the bed file containing 
                      CDSs annotation (chr, start, end).')
  parser$add_argument('outPath', metavar = 'output', type = "character", 
                      nargs = 1, help = 'Full path to the output folder.')
  parser$add_argument('path2OCFiles', metavar = 'OCPFiles', type = "character", 
                      nargs = 1, help = 'Full path to the folder where 
                      atac-seq/dnase-seq peak files are.')
  parser$add_argument('--path2K4me1Files', dest = 'path2K4me1Files', 
                      action = 'store', default = NULL, 
                      help = 'Full path to the folder where H3K4me1 ChiP-seq 
                      peaks files are.')
  parser$add_argument('--path2K27acFiles', dest = 'path2K27acFiles', 
                      action = 'store', default = NULL, 
                      help = 'Full path to the folder where H3K27ac ChiP-seq 
                      peaks files are.')
  parser$add_argument('--path2K4me3Files', dest = 'path2K4me3Files', 
                      action = 'store', default = NULL, 
                      help = 'Full path to the folder where H3K4me3 ChIP-seq 
                      peaks files are.')
  parser$add_argument('--oldAnnotation', dest = 'oldAnnotation', 
                      action = 'store', default = NULL, help = 'Full path to 
                      the gff3 file with previous regulatory annotation.')
  parser$add_argument('--ocPattern', dest = 'ocPattern', action = 'store',
                      default = '', help = 'Pattern characterizing the subset of 
                      atac-seq/dnase-seq bed files to be used.')
  parser$add_argument('--k4me1Pattern', dest = 'k4me1Pattern', 
                      action = 'store', default = '', 
                      help = 'Pattern characterizing the subset of H3K4me1 
                      bed files to be used.')
  parser$add_argument('--k27acPattern', dest = 'k27acPattern',
                      action = 'store', default = '',
                      help = 'Pattern characterizing the subset of 
                      H3K27ac bed files to be used.')
  parser$add_argument('--k4me3Pattern', dest = 'k4me3Pattern',
                      action = 'store', default = '',
                      help = 'Pattern characterizing the subset of 
                      H3K4me3 bed files to be used.')
  parser$add_argument('--TSSwindowUp', dest = 'TSSwindowUp', action = 'store', 
                      default = 90, help = 'Window size (bp) upstream TSSs for 
                      promoters identification (default 90).')
  parser$add_argument('--TSSwindowDown', dest = 'TSSwindowDown', 
                      action = 'store', default = 10, help = 'Window size (bp)
                      downstream TSSs for promoters identification 
                      (default 10).')
  parser$add_argument('--coreWUp', dest = 'coreWUp', action = 'store', 
                      default = 490, help = 'Upstream core size (bp) for 
                      potential promoters (default 490).')
  parser$add_argument('--coreWDown', dest = 'coreWDown', action = 'store', 
                      default = 10, help = 'Downstream core size (bp) for
                      potential promoters (default 10).')
  parser$add_argument('--featureOverlap', dest = 'featureOverlap', 
                      action = 'store', default = 0.2, help = 'Fraction of 
                      overlapping between RFs and open chromatin peaks
                      for activity annotation (default 0.2).')
  parser$add_argument('--histoneOverlap', dest = 'histoneOverlap', 
                      action = 'store', default = 0.5, help = 'Fraction of 
                      overlap between histone and open chromatin peaks, for 
                      features identification (default 0.5).')
  parser$add_argument('--OCRsHistoneOverlap', dest = 'OCRsHistoneOverlap', 
                      action = 'store', default = 0.5, help = 'Fraction of 
                      overlap between cOCRs and histone peaks for enhancers 
                      identification (default 0.5).')
  parser$add_argument('--cOCRsCDSOverlap', dest = 'cOCRsCDSOverlap', 
                      action = 'store', default = 0.1, help = 'Fraction of 
                      overlap between cOCRs and CDSs for enhancers 
                      identification (default 0.1).')
  parser$add_argument('--maxDistOCR', dest = 'maxDistOCR', action = 'store', 
                      default = 100, help = 'Maximum distance between 
                      peaks to be merged when defining cOCRs (default 100).')
  parser$add_argument('--maxDistH', dest = 'maxDistH', action = 'store', 
                      default = 100, help = 'Maximum distance between histone 
                      peaks to be merged (default 100).')  
  parser$add_argument('--maxPromExtLen', dest = 'maxPromExtLen', 
                      action = 'store', default = 1000, help = 'Maximum size 
                      for promoters extended regions (default 1000).')
  parser$add_argument('--minOCRLen', dest = 'minOCRLen', action = 'store',
                      default = 100, 
                      help = 'Minimum size for cOCRs (default 100).')
  parser$add_argument('--maxEnhancerLen', dest = 'maxEnhancerLen',
                      action = 'store', default = 2500, help = 'Maximum length 
                      for enhancers (default 2500).')  
  parser$add_argument('--maxOCRLen', dest = 'maxOCRLen', action = 'store', 
                      default = 2500, help = 'Maximum length for open
                      chromatin regions (default 2500).')
  parser$add_argument('--exportActivity', dest = 'exportActivity', 
                      action = 'store', default = FALSE, help = 'Should 
                      activity files be generated?(default FALSE)')
  parser$add_argument('--path2EpigenomesOCFiles', dest = 'EpiOCFiles', 
                      action = "store", default = '',
                      help = 'Full path to the folder containing the open 
                      chromatin peak files to be used for activity annotation.')
  parser$add_argument('--epiPattern', dest = 'epiPattern',
                      action = 'store', default = '',
                      help = 'Pattern characterizing the subset of 
                      bed files to be used for activity annotation.')
  
  return(parser)
}

loadPeaks <- function(path2Files, bedFiles, chrs){
  allPeaks <- list()
  # Read bed files
  for(bf in bedFiles){
    peaks <- tryCatch(read.delim(paste(path2Files, bf, sep = "/"), 
                                 header = F),
                      warning = function(w){
                        print(paste(bf, 'does not contain any peak', 
                                    sep = " "))
                        return(NULL)}
    )
    # Keep only peaks in chromosomes of interest
    peaks %>%
      mutate(V1 = as.character(V1)) %>%
      filter(V1 %in% chrs) -> allPeaks[[bf]]
  }
  return(allPeaks)
}

loadEpiPeaks <- function(path2EpigenomesOCFiles, epiBedFiles, chrs, assembly){
  # Epigenome name is embedded in the file names
  epigenomes <- do.call(rbind, strsplit(epiBedFiles, split = assembly))[,2]
  epigenomes <- do.call(rbind, strsplit(epigenomes, split = '-seq-'))[,2]
  # Each epigenome could have more than one linked open-chromatin peak file
  epigenomes <- unique(do.call(rbind, strsplit(epigenomes, split = '_'))[,1])
  epiPeaks <- list()
  for(epi in epigenomes){
    bedFiles <- epiBedFiles[grepl(paste(epi, '_', sep =''), epiBedFiles)]
    peaks <- NULL
    # Load all the peak files and merge overlapping peaks for having a set of
    # unique open chromatin regions 
    for(i in 1:length(bedFiles)){
      peaks %>%
        bind_rows(read.delim(paste(path2EpigenomesOCFiles, bedFiles[i],
                                            sep = "/"), header = F) %>%
                    mutate(V1 = as.character(V1))) %>%
        bt.sort() %>%
        bt.merge() %>%
        filter(V1 %in% chrs) -> peaks
    }
    epiPeaks[[epi]] <- peaks
  }
  return(epiPeaks)
}

defineUnionPeaks <- function(peaksList, d = NULL){
  # Merge nearby peaks up to d bp
  bind_rows(peaksList) %>%
    dplyr::select(V1, V2, V3) %>%
    bt.sort() %>%
    bt.merge(d = d) -> uPeaks
  return(uPeaks)
}

addingTSSsInfo <- function(annotatedTSSsAndInfo, TSSwindowUp, TSSwindowDown,
                           coreWUp, coreWDown){
  # Create TSSs windows and core region useful for promoter identification
  annotatedTSSsAndInfo %>% 
    select(V1, V2, V3, V4) %>%
    unique() -> annotatedTSSs
  names(annotatedTSSs) <- c('TSSChr', 'TSSStart', 'TSSEnd', 'TSSStrand')
  # The strand need to be taken into account
  annotatedTSSs %>%
    mutate(wStart = if_else(TSSStrand == '+',
                            TSSStart - TSSwindowUp, TSSStart - TSSwindowDown ),
           wEnd= if_else(TSSStrand == '+',
                         TSSEnd + TSSwindowDown, 
                         TSSEnd + TSSwindowUp )) -> annotatedTSSs
  annotatedTSSs %>%
    mutate(cStart = if_else(TSSStrand == '+',
                            TSSStart - coreWUp, TSSStart - coreWDown ),
           cEnd = if_else(TSSStrand == '+',
                          TSSEnd + coreWDown, 
                          TSSEnd + coreWUp )) -> annotatedTSSs
  annotatedTSSs %>%
    mutate(across(c(wStart, cStart), 
                  function(x){if_else( x < 0, 0, x )})) -> annotatedTSSs
  
  annotatedTSSs %>%
    dplyr::select(TSSChr, wStart, wEnd, TSSStrand, TSSStart, TSSEnd, 
                  cStart, cEnd) -> annotatedTSSs
  return(annotatedTSSs)
}
definingUPsInTSS <- function(UPs, annotatedTSSs){
  # Identify the open chromatin UPs overlapping TSSs windows
  UPs %>%
    bt.intersect(annotatedTSSs, wa = T, wb = T) -> UPsTSS
  
  names(UPsTSS) <- c(c('upChr', 'upStart', 'upEnd'), names(annotatedTSSs))
  # Trimming/extending UPs-TSSs for excluding regions up TSS wEnd from
  # the TSS defining final coordinates
  UPsTSS %>% 
    mutate(upNewStart = if_else(TSSStrand == '+', 
                                pmin(upStart,  wStart), 
                                pmin(pmax(upStart, wStart), wStart)),
           upNewEnd =  if_else(TSSStrand == '+',
                               pmin(pmax(upEnd, wEnd), wEnd),
                               pmax(upEnd, wEnd))) -> UPsTSS
  
  # updating core region to the experimental information
  UPsTSS %>%
    mutate(cStart = if_else(TSSStrand == '+' & (cStart < upNewStart),
                            upNewStart, cStart),
           cEnd = if_else(TSSStrand == '-' & (cEnd > upNewEnd),
                          upNewEnd, cEnd)) -> UPsTSS
  
  return(UPsTSS)
}
promotersIdentification <- function(UPsTSS, TSSwindowDown, maxPromExtLen){
  UPsTSS %>%
    dplyr::select(upChr, upStart, upEnd) %>%
    unique() -> uniqueUPsTSS
  
  promoters <- data.frame(chr = NULL, core_start = NULL, core_end=NULL,
                          extended_start = NULL, extended_end = NULL)
  cOCRs <- data.frame(chr = NULL, start = NULL, end=NULL)
  # each UP overlapping at least a TSS is evaluated
  for(i in 1:nrow(uniqueUPsTSS)){
    up <- uniqueUPsTSS[i,]
    UPsTSS %>%
      filter(upChr == up[ ,'upChr'] &
               upStart == up[ ,'upStart'] &
               upEnd == up[ ,'upEnd']) %>%
      arrange(TSSChr, TSSStart, TSSEnd) -> upTss
    mergedCores <- NULL
    if(nrow(upTss) > 1){
      r1 <- upTss[1 , c('TSSChr', 'cStart', 'cEnd', 'TSSStrand')]
      for(tss in 2:nrow(upTss) ){
        r2 <- upTss[tss , c('TSSChr', 'cStart', 'cEnd','TSSStrand')]
        
        if(r1[1, 'TSSStrand'] == '-' & r2[1, 'TSSStrand'] == '+' & 
           r1[1, 'cEnd'] >= r2[1, 'cEnd'] & 
           r1[1, 'cStart'] >= r2[1, 'cStart'] ){
          mc <- bt.intersect(a = r1[,1:3], b = r2[,1:3])
        }else{
          if(r1[1, 'TSSStrand'] == '-' & r2[1, 'TSSStrand'] == '+' & 
             r1[1, 'cEnd'] > r2[1, 'cEnd'] & 
             r2[1, 'cStart'] > r1[1, 'cStart']){
            #correction for mergedCore previously extended for a '-' TSS
            r1[1, 'cEnd'] <- r2[1, 'cEnd']
          }
          if(r2[1, 'cStart'] < r1[1, 'cStart']){
            #correction for mergedCore previously trimmed for a '-' TSS
            r2[1, 'cStart'] <- r1[1, 'cStart']
          }
          mc <- bt.merge(bt.sort(rbind(r1[, 1:3], r2[, 1:3])))
        }
        names(mc) <- c('TSSChr', 'cStart','cEnd')
        if(nrow(mc) > 1){
          if(tss == 2){
            mergedCores <- rbind(mergedCores, mc[1, ])
          }
          mergedCores <- rbind(mergedCores, mc[nrow(mc), ])
        }else{
          if(tss == 2){
            mergedCores <- rbind(mergedCores, mc[1, ])
          }else{
            mergedCores[nrow(mergedCores), ] <- mc 
          }
        }
        r1 <- cbind(mc[nrow(mc), ], TSSStrand = r2[1, 'TSSStrand'])
      }
    }else{
      mergedCores <- upTss[1, c('TSSChr', 'cStart', 'cEnd')]
    }
    if(nrow(mergedCores) == 1 ){# up will define a single promoter
      proms <-data.frame(chr = mergedCores[, 'TSSChr'], 
                         core_start = as.numeric(upTss[1, 'cStart']), 
                         core_end = as.numeric(upTss[nrow(upTss), 'cEnd']),
                         extended_start = as.numeric(upTss[1, 'upNewStart']), 
                         extended_end = as.numeric(upTss[nrow(upTss), 
                                                         'upNewEnd']))
    }else{# more than one promoter core
      proms <- data.frame(chr = NULL, 
                          core_start = NULL, 
                          core_end = NULL,
                          extended_start = NULL, 
                          extended_end = NULL)
      for(j in 1:nrow(mergedCores)){
        case1 <- upTss[, 'cStart' ] >= mergedCores[j, 'cStart'] & 
          upTss[, 'cEnd' ] <= mergedCores[j, 'cEnd' ]
        case2 <- upTss[, 'cStart' ] == mergedCores[j, 'cStart'] & 
          upTss[, 'cEnd' ] >= mergedCores[j, 'cEnd' ]
        case3 <- upTss[, 'cStart' ] <= mergedCores[j, 'cStart'] & 
          upTss[, 'cEnd' ] == mergedCores[j, 'cEnd' ]
        cond1 <- upTss[, 'TSSStart'] >= mergedCores[j, 'cStart'] & 
          upTss[, 'TSSEnd' ] <= mergedCores[j, 'cEnd' ]
        
        overlapped <- upTss[cond1 & (case1 | case2 | case3), ]
        proms <- rbind(proms, 
                       data.frame(chr = mergedCores[j, 'TSSChr'], 
                                  core_start = as.numeric(mergedCores[j, 'cStart']), 
                                  core_end = as.numeric(mergedCores[j, 'cEnd' ]),
                                  extended_start = as.numeric(overlapped[1, 'upNewStart']), 
                                  extended_end = as.numeric(overlapped[nrow(overlapped),
                                                                       'upNewEnd'])))
      }
      proms <- as.data.frame(proms)
      for(p in 1:(nrow(proms) -1)){
        if(proms[p, 'extended_end'] == up[, 'upEnd'] & 
           proms[p + 1, 'extended_start'] == up[, 'upStart'] ){
          # correct extended regions of consecutive promoters
          proms[p, 'extended_end'] <- proms[p, 'core_end'] + 
            round((proms[p + 1, 'core_start'] - proms[p, 'core_end']) / 2)
          proms[p + 1, 'extended_start'] <- proms[p, 'extended_end'] 
        }
        if(proms[p, 'extended_start'] == proms[p + 1, 'extended_start']){
          # Two consecutive Tss in '+', prom(p) should end before next start
          proms[p + 1, 'extended_start'] <- proms[p, 'extended_end'] 
        }
        if(proms[p, 'extended_end'] >= proms[p + 1, 'extended_start']){
          if(proms[p, 'extended_start'] < proms[p + 1, 'extended_start'] & 
             proms[p, 'extended_end'] > proms[p, 'core_end']){
            # Two consecutive Tss in '-', prom(p) should be trim before next start
            proms[p, 'extended_end'] <- proms[p + 1, 'extended_start']
          }else{
            # More than two consecutive Tss in '+',prom(p +1) should start after prom(p)
            proms[p + 1 , 'extended_start'] <- proms[p, 'extended_end']
          }
        }
        if(proms[p + 1, 'extended_start'] - proms[p, 'extended_end'] > 1){
          #First promoter in + next in -: intermediate ocr
          cOCRs <- rbind(cOCRs, data.frame(chr = proms[p,'chr'],
                                           start = proms[p, 'extended_end'] ,
                                           end = proms[p + 1, 'extended_start']))
        }
      }
    }
    proms <- as.data.frame(proms)
    promsToAdd <- proms
    if(i > 1){
      #checking potential overlapping caused by the extension of UPs overlapping
      # TSSs windows
      lastP <- promoters[nrow(promoters), ]
      if((lastP[, 'chr'] == proms[1, 'chr']) & 
         (lastP[, 'extended_end'] > proms[1, 'extended_start'])){
        if(lastP[, 'core_end'] > proms[1, 'core_start']){
          # cores overlapping, keep only one promoter merging overlapping cores
          promsStrand <- upTss[upTss[, 'cEnd'] == proms[1, 'core_end'], ]
          promsStrand <- promsStrand[nrow(promsStrand), 'TSSStrand']
          if(lastTSSstrand == promsStrand){
            lastP[, 'core_end'] <- proms[1, 'core_end']
            lastP[, 'extended_end'] <- proms[1, 'extended_end']
          }else{
            if(lastTSSstrand == '-'){
              # it can occur because of the extension there is a '+' TSS
              # that didn't overlap with the UP but that is in the extended 
              # region. In this case the promoter should not be extended and the 
              # UP will be considered for defining OCRs
              oldBS <- lastP[, 'extended_start']
              oldBE <- lastP[, 'extended_end']
              lastP[, 'core_end'] <- proms[1, 'core_end']
              lastP[, 'extended_end'] <- proms[1, 'extended_end']
              # because we maybe trim the lastP
              if(oldBS < lastP[, 'extended_start']){
                cOCRs <- rbind(cOCRs, data.frame(chr = lastP[, 'chr'],
                                                 start = oldBS,
                                                 end = lastP[, 'extended_start']))
              }
              if(oldBE > lastP[, 'extended_end']){
                cOCRs <- rbind(cOCRs, data.frame(chr = lastP[, 'chr'],
                                                 start = lastP[, 'extended_end'],
                                                 end = oldBE))
              }
            }else{
              if((lastP[, 'core_start'] < proms[1, 'core_start']) & 
                 (proms[1, 'core_start'] < 
                  (lastP[, 'core_end'] - 2*TSSwindowDown))){
                #add OCR upstream
                cOCRs <- rbind(cOCRs, data.frame(chr = lastP[, 'chr'],
                                                 start = lastP[, 'extended_start'],
                                                 end = proms[1, 'core_start']))  
                lastP[, 'core_start'] <- proms[1, 'core_start']
                lastP[, 'extended_start'] <- proms[1, 'extended_start']
              }
              # check if there are some TSSs downstream the lastP, we need
              # to redefine the prom
              upTss %>%
                filter(cStart > (lastP[, 'core_end'] - 2*TSSwindowDown) &
                         cEnd <= proms[1, 'core_end']) -> extupTss
              if((nrow(extupTss) > 0)){
                if((min(extupTss[, 'cStart']) > lastP[, 'core_end'])){
                  proms <- rbind(proms, proms[1,])  
                  proms[nrow(proms), 'core_start'] <- min(extupTss[, 'cStart'])
                  proms[nrow(proms), 'extended_start'] <- min(extupTss[, 'cStart'])
                  if(extupTss[1, 'upStart'] < min(extupTss[, 'cStart'])){
                    cOCRs <- rbind(cOCRs, 
                                   data.frame(chr = extupTss[1, 'upChr'],
                                              start = max(lastP[, 'extended_end'],
                                                          extupTss[1, 'upStart']),
                                              end = min(extupTss[, 'cStart'])))
                    
                  }
                  proms <- proms[order(proms[, 'core_start']),]
                }else{
                  # extend the current lastP
                  lastP[, 'core_end'] <- proms[1, 'core_end']
                  lastP[, 'extended_end'] <- proms[1, 'extended_end']
                }
              }else{
                #to ensure last strand is the '+'
                upTss[upTss[, 'cEnd'] == proms[1, 'core_end'], 
                      'TSSStrand'] <- lastTSSstrand
              }
            }
          }
          promsToAdd <- proms[-1,]
          proms[1,] <- lastP
        }else{ 
          # fixing overlapping of extended regions 
          lastP[, 'extended_end'] <- proms[1, 'extended_start'] + 
            round((proms[1, 'extended_start'] - 
                     lastP[, 'extended_end']) / 2)
          proms[1, 'extended_start'] <- lastP[, 'extended_end'] 
          promsToAdd <- proms
        }
        promoters[nrow(promoters), ] <- lastP
        # remove internal cOCRs that are now overlapping extended promoters
        cocr2remove <- bt.intersect(cOCRs, rbind(lastP, promsToAdd), f = 1)
        if(nrow(cocr2remove) >=1){
          for(ci in 1:nrow(cocr2remove)){
            cOCRs %>%
              filter((chr != cocr2remove[ci, 1]) |
                       (start != cocr2remove[ci, 2] & 
                          end != cocr2remove[ci, 3])) -> cOCRs
          }
        }
      }   
    }
    lastTSSstrand <- upTss[nrow(upTss), 'TSSStrand']
    promoters <- rbind(promoters, promsToAdd)
    if(nrow(proms) >= 1){
      if(proms[1, 'extended_start'] > up[, 'upStart']){ # OCR upstream
        cOCRs <- rbind(cOCRs, data.frame(chr = up[, 'upChr'],
                                         start = up[, 'upStart'],
                                         end = proms[1, 'extended_start']))
      }
      if(proms[nrow(proms), 'extended_end'] < up[, 'upEnd']){ # OCR downstream
        cOCRs <- rbind(cOCRs, data.frame(chr = up[, 'upChr'],
                                         start = proms[nrow(proms), 'extended_end']  ,
                                         end = up[, 'upEnd']))
      }
      
    }
  }
  
  #because some short UP could be extended because overlapping a '+' TSS but 
  # then they were identified as OCR
  uniqueUPsTSS %>% 
    bt.intersect(b = cOCRs, f = 1) -> upsIncOCRs
  cOCRs %>% 
    bt.intersect(b = upsIncOCRs, v = T) -> cOCRs
  names(cOCRs) <- c('chr', 'start', 'end')
  # If it was specified, trim promoter bounds to maxPromExtLen
  if(!is.null(maxPromExtLen)){
    promoters %>% 
      filter((core_start - extended_start) > maxPromExtLen) %>%
      dplyr::select(chr, extended_start, core_start) %>%
      mutate(core_start = core_start - maxPromExtLen) %>%
      rename('start' = extended_start, 'end' = core_start) %>%
      bind_rows(promoters %>% 
                  filter((extended_end - core_end) > maxPromExtLen) %>%
                  dplyr::select(chr, core_end, extended_end) %>%
                  mutate(core_end = core_end + maxPromExtLen) %>%
                  rename('start' = core_end, 'end' = extended_end)) %>%
      mutate(chr = as.character(chr)) %>%
      bind_rows(cOCRs %>% mutate(chr = as.character(chr))) -> cOCRs
    promoters %>% 
      mutate(extended_start = if_else((core_start - extended_start) > maxPromExtLen,
                                      core_start - maxPromExtLen, 
                                      extended_start),
             extended_end = if_else((extended_end - core_end) > maxPromExtLen,
                                    core_end + maxPromExtLen, 
                                    extended_end)) -> promoters
  }
  return(list('promoters'=promoters, 'cOCRs'=cOCRs, 'upsIncOCRs' = upsIncOCRs))
}

OCRsIdentification <- function(UPsnoTSS, cOCRs, promoters, maxDistOCR, 
                               minOCRLen){
  # combine all the cOCRs and remaining UPs fragments
  cOCRs %>% bt.sort() %>% bt.merge() %>%
    mutate(V1 = as.character(V1)) %>%  
    rename( 'chr' = V1, 'start' = V2, 'end' = V3) %>% 
    bind_rows(UPsnoTSS) %>% bt.sort() -> OCRs
  # This step is aimed to merge OCRs that are up to maxDistOCR checking
  # that no overlapping with promoters are produced
  mcOCRs <- bt.merge(OCRs, d = maxDistOCR) # merge OCRs,
  imOCRs <- which(bt.coverage(mcOCRs, OCRs)[,4] > 1) # mOCRs index
  # for each merged regions, we should check no promoters are overlapping
  imToExclude <- NULL
  ocrToAdd <- NULL
  for(i in imOCRs){
    mr <- mcOCRs[i, ]
    if(any(promoters[, 'chr'] == mr[,1] & 
           promoters[, 'extended_start'] >= mr[, 2] &
           promoters[, 'extended_end'] <= mr[, 3])){
      pI <- promoters[promoters[, 'chr'] == mr[, 1] & 
                        promoters[, 'extended_start'] >= mr[, 2] &
                        promoters[, 'extended_end'] <= mr[, 3],]
      # if there is a promoter in the middle, the OCRs shouldn't be merged
      imToExclude <- c(imToExclude, i)
      # and the OCRs regions should be added to the merged OCRs
      ocrs <- OCRs[OCRs[, 1] == mr[, 1] & 
                     OCRs[, 2] >= mr[, 2] & 
                     OCRs[, 3] <= mr[, 3], ]
      #because some of the OCRs can still be merged
      ocrsToM <- (ocrs[, 1] == unique(pI[, 'chr']) & 
                    ocrs[, 3] < min(pI[, 'extended_start']) - maxDistOCR) | 
        (ocrs[, 1] == unique(pI[, 'chr']) & 
           ocrs[, 2] > max(pI[, 'extended_end']))
      ocrsM <- bt.merge(ocrs[ocrsToM, ], d = maxDistOCR)
      ocrs <- rbind(ocrs[!ocrsToM, ], ocrsM)
      ocrToAdd <- rbind(ocrToAdd, ocrs)
    }
  }
  mcOCRs <- mcOCRs[!rownames(mcOCRs) %in% imToExclude,]
  mcOCRs <- rbind(mcOCRs, ocrToAdd)
  mcOCRs <- bt.sort(mcOCRs)
  # just for compatibility of data.frames, OCRs will have core and extended
  # regions but both have same length
  mcOCRs <- cbind(mcOCRs, mcOCRs[,2],mcOCRs[,3])
  colnames(mcOCRs) <- colnames(promoters)
  # Filter artifact OCRs
  if(!is.null(minOCRLen)){
    mcOCRs <- mcOCRs[mcOCRs[ ,'extended_end'] - 
                       mcOCRs[, 'extended_start'] >= minOCRLen, ] 
  }
  return(mcOCRs)
}
enhancersIdentification <- function(mcOCRs, minOCRLen,
                                    H3K4me1OCUPs, H3K27acOCUPs, 
                                    annotatedCDS, cOCRsCDSOverlap,
                                    OCRsHistoneOverlap, maxEnhancerLen){
  # If present, filter mcOCRs by overlapping with H3K4me1OCUPs to decide which 
  # ones could be defined as enhancers
  if(!is.null(H3K4me1OCUPs)){
    mcOCRs %>%
      dplyr::select(chr, extended_start, extended_end) %>%
      bt.intersect(b = H3K4me1OCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T,
                   u = T) -> enhancers
  }else{
    enhancers <- NULL
  }
  # If present, filter mcOCRs by overlapping with H3K27acOCUPs to decide which 
  # ones could be defined as enhancers
  if(!is.null(H3K27acOCUPs)){
    mcOCRs %>%
      dplyr::select(chr, extended_start, extended_end) %>%
      bt.intersect(b = H3K27acOCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T,
                   u = T) -> enhancers27ac
    # combine previously identified enhancers
    enhancers %>%
      bind_rows(enhancers27ac) %>%
      unique() -> enhancers
  }
  # Maximum length for enhancers
  if(!is.null(maxEnhancerLen)){
    shrinkedE <- NULL
    # Try first to shrink the long enhancers to H3K27ac peaks when possible
    if(!is.null(H3K27acOCUPs)){
      enhancers %>%
        filter((V3 - V2) > maxEnhancerLen) %>%
        bt.intersect(b = H3K27acOCUPs, e = T,
                     f = OCRsHistoneOverlap, 
                     F = OCRsHistoneOverlap) -> shrinkedE
    }
    enhancers %>%
      bind_rows(shrinkedE) %>%
      filter((V3 - V2) <= maxEnhancerLen) -> enhancers
  }
  # discard enhancers overlapping CDS
  enhancers %>%
    bt.intersect(b = annotatedCDS, wao = TRUE) %>%
    group_by(V1, V2, V3) %>%
    summarise(cumOV = sum(V7), .groups = 'drop') %>%
    mutate(maxOV = (V3 - V2) * cOCRsCDSOverlap) %>%
    filter(cumOV < maxOV) %>%
    dplyr::select(V1, V2, V3) %>%
    as.data.frame() -> enhancers 
  enhancers <- cbind(enhancers, enhancers[,2], enhancers[,3])
  colnames(enhancers) <- colnames(mcOCRs)
  
  # remove enhancers shorter than minOCRLen
  if(!is.null(minOCRLen)){
    enhancers %>%
      filter(extended_end - extended_start >= minOCRLen) -> enhancers
  }
  # Refine final set of OCRs by filtering mcOCRs that overlap enhancers or
  # their histone marks
  mcOCRs %>%
    bt.intersect(enhancers, wa = TRUE, v = TRUE) -> mcOCRs
  
  if(!is.null(H3K4me1OCUPs)){
    mcOCRs %>%
      bt.intersect(b = H3K4me1OCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T, 
                   v = T) -> mcOCRs
  }
  
  if(!is.null(H3K27acOCUPs)){
    mcOCRs %>%
      bt.intersect(b = H3K27acOCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T, 
                   v = T) -> mcOCRs
  }
  colnames(mcOCRs) <- colnames(enhancers)
  return(list(enhancers = enhancers, mcOCRs = mcOCRs))
}

emarsIdentification <- function(ocUPs, H3K4me1OCUPs, H3K27acOCUPs,
                                H3K4me3OCUPs, OCRsHistoneOverlap, minOCRLen){
  # EMARs are open chromatin regions with histone modifications
  if(!is.null(H3K4me1OCUPs)){
    ocUPs %>% 
      bt.intersect(b = H3K4me1OCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T, 
                   u = T) -> emars
  }else{
    emars <- NULL
  }
  if(!is.null(H3K27acOCUPs)){
    ocUPs %>% 
      bt.intersect(b = H3K27acOCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T, 
                   u = T) -> emars27ac
    emars %>%
        bind_rows(emars27ac) %>%
        unique() -> emars
  }
  if(!is.null(H3K4me3OCUPs)){
    ocUPs %>% 
      bt.intersect(b = H3K4me3OCUPs, e = T,
                   f = OCRsHistoneOverlap, F= OCRsHistoneOverlap, wa = T, 
                   u = T) -> emars4me3
    emars %>%
      bind_rows(emars4me3) %>%
      unique() -> emars
  }
  emars %>%
    bt.sort() %>%
    bt.merge() -> emars
  # Remove EMARs shorter than minOCRLen
  if(!is.null(minOCRLen)){
    emars %>%
      filter( V3 - V2 >= minOCRLen) -> emars
  }
  emars <- cbind(emars, emars[,2], emars[,3])
  names(emars) <- c('chr', 'core_start', 'core_end', 'extended_start', 
                    'extended_end')
  return(emars)
}
createRFsGFF <- function(rfs, chrs, ESpePre, startID, previousRegBuild, 
                         annotatedTSSsAndInfo, 
                         featuresToMerge = c('promoter', 'enhancer', 
                                             'open_chromatin_region', 'EMAR'),
                         source = 'Ensembl'){
  # export features in gff3 format
  rfs %>%
    mutate(seqid = factor(chr, levels = chrs, ordered = TRUE),
           source = source,
           type = type,
           start = core_start + 1,
           end = core_end,
           score = '.',
           strand = '.', 
           phase = '.',
           extended_start = extended_start + 1) %>%
    dplyr::select(seqid, source, type, start, end, score, strand, phase, 
           extended_start, extended_end) %>%
    arrange(type, seqid, start) -> rfsGFF
  
  rownames(rfsGFF) <- paste(rfsGFF[,'type'], rfsGFF[,'seqid'],
                            paste(rfsGFF[,'extended_start'], 
                                  rfsGFF[,'extended_end'],
                                  sep = '-'), sep = ":")
  ID <- rep('', nrow(rfsGFF))
  names(ID) <- rownames(rfsGFF) 
  if(!is.null(previousRegBuild)){
    #Identify old features mapping one to one with new features to inherit the 
    # ID
    # Minimum overlapping and ignoring matching between feature types
    # If TFBS/CTCF are annotated, then consider filter to overlap just
    # c( 'enhancer','open_chromatin_region','promoter')
    # EMARs are analyzed separately too
    
    previousRegBuild %>%
      filter(type %in% featuresToMerge) -> previousRegBuild2Merge
    if('EMAR' %in% featuresToMerge){
      previousRegBuild2Merge %>%
        filter(type == 'EMAR') -> previousEMARs
      previousRegBuild2Merge %>%
        filter(type != 'EMAR') -> previousRegBuild2Merge
      
      featuresToMerge <- featuresToMerge[featuresToMerge != 'EMAR']
      rfsGFF %>%
        filter(type == 'EMAR') %>%
        dplyr::select(seqid, extended_start, extended_end, type) %>%
        bt.intersect(previousEMARs, wa = TRUE, wb = TRUE) %>%
        mutate(loc = paste(V4, V1, paste(V2, V3, sep ='-'), 
                           sep = ':')) -> EMARsOne2One 
      # Remove the one-to-many cases
      EMARsOne2One %>% 
        group_by(V9) %>% 
        summarise(n=n()) %>% 
        filter(n > 1) %>%
        as.data.frame() -> EMARsOne2Many
      # Remove the many-to-one cases
      EMARsOne2One %>% 
        group_by(loc) %>% 
        summarise(n=n()) %>% 
        as.data.frame() %>% 
        filter( n > 1)-> EMARsMany2One
      EMARsOne2One %>% 
        filter(!(V9 %in% EMARsOne2Many[,1]) &
                 !(loc %in% EMARsMany2One[,1]) )  -> EMARsOne2One
      ID[EMARsOne2One[,'loc']] <- EMARsOne2One[, 9]
    }
    
    rfsGFF %>%
      filter(type %in% featuresToMerge) %>%
      dplyr::select(seqid, extended_start, extended_end, type) %>%
      bt.intersect(previousRegBuild2Merge, wa = TRUE, wb = TRUE) %>%
      mutate(loc = paste(V4, V1, paste(V2, V3, sep ='-'), 
                         sep = ':')) -> rfsGFFOne2One 
    # Remove the one-to-many cases
    rfsGFFOne2One %>% 
      group_by(V9) %>% 
      summarise(n=n()) %>% 
      filter(n > 1) %>%
      as.data.frame() -> oneToMany
    
    # Remove the many-to-one cases
    
    rfsGFFOne2One %>% 
      group_by(loc) %>% 
      summarise(n=n()) %>% 
      filter(n > 1) %>%
      as.data.frame() -> manyToOne
    
    rfsGFFOne2One %>% 
      filter(!(V9 %in% oneToMany[,1]) &
               !(loc %in% manyToOne[,1]) )  -> rfsGFFOne2One
    ID[rfsGFFOne2One[,'loc']] <- rfsGFFOne2One[, 9] 
  }
  nNewFeat <- length(which(ID == ''))
  IDNewFeat <- 1:nNewFeat + startID - 1
  IDNewFeat <- sprintf(paste(ESpePre, "R%011d", sep = ""), IDNewFeat) 
  ID[ID == ''] <- IDNewFeat
  ID <- paste('ID=', ID, sep= '')
  # we only report extended regions for promoters
  BS <- if_else(rfsGFF[,'type'] == 'promoter',
                 paste(';extended_start=', rfsGFF[,'extended_start'], sep= ''),
                 '')
  BE <- if_else(rfsGFF[,'type'] == 'promoter',
                paste(';extended_end=', rfsGFF[,'extended_end'], sep= ''),
                '')
  # adding gene info for promoters
  rfs %>% 
    filter(type == 'promoter') %>%
    dplyr::select(chr, extended_start, extended_end) %>%
    bt.intersect(annotatedTSSsAndInfo,
                 wa = TRUE, wb = TRUE) -> rfsProm
  if(ncol(rfsProm) == 10){
    rfsProm %>% mutate(V10 = if_else(V10 == V8, '', V10)) -> rfsProm
  }else{
    rfsProm %>%
      bind_cols('V10' = '') -> rfsProm
  }
  rfsProm %>%  
    dplyr::select(-c(V4, V5, V6, V7)) %>%
    unique() %>%
    group_by(V1, V2, V3) %>%
    summarise(gene_id = paste(V8, collapse =','),
              gene_name = paste(V10, collapse =','),
              gene_biotype = paste(V9, collapse =','), .groups = 'drop') %>%
    mutate(GI = if_else(str_remove_all(gene_name,',') == '',
                        paste(';gene_id=', gene_id, sep = ''),
                        paste(paste(';gene_id=', gene_id, sep = ''),
                              paste('gene_name=', gene_name, sep = ''), 
                              sep = ';'))) %>%
    mutate(GI = paste(GI, paste('gene_biotype=', gene_biotype, sep = ''), 
                      sep = ';')) %>%
    dplyr::select(V1, V2, V3, GI) %>%
    as.data.frame() -> GIInfo
             
  rownames(GIInfo) <- paste('promoter', GIInfo[,1], 
                            paste(GIInfo[,2] + 1,
                                  GIInfo[,3], sep = '-'), sep =':')
  GI <- rep('', nrow(rfsGFF))
  names(GI) <- rownames(rfsGFF)
  GIInfo %>% pull(GI) -> GI[rownames(GIInfo)]
  
  rfsColors <- c('promoter' = '#ff0000',
                 'enhancer' = '#faca00',
                 'open_chromatin_region' = '#d9d9d9',
                 'CTCF_binding_site' = '#40e0d0',
                 'EMAR' = '#004d40')
  FC <- rfsColors[rfsGFF[,'type']]
  FC <- paste(';color=', FC, sep= '')
  rfsGFF[, 'attributes'] <- paste(ID, BS, BE, GI, FC, sep= '')
  rfsGFF <- rfsGFF[, !(colnames(rfsGFF) %in% c('extended_start','extended_end'))]
  rownames(rfsGFF) <- paste(rfsGFF[,'type'], rfsGFF[,'seqid'],
                            paste(rfsGFF[,'start'] - 1, 
                                  rfsGFF[,'end'],
                                  sep = '-'), sep = ":")
  rfsGFF %>%
    mutate(seqid = factor(seqid, levels = chrs, ordered = TRUE)) %>%
    arrange(seqid, start) -> rfsGFF
  return(rfsGFF)
}
rfsInEpigenomes <- function(rfs, epiPeaks, annotatedTSSs, featureOverlap, eVal){
  # The overlapping is computed differently for promoters than for enhancers and
  # open chromatin regions. For the latter, we request minimum 
  # overlapping fraction (featureOverlap)
  
  # exclude EMARs
  rfs %>%
    filter(type != 'EMAR') -> rfs
  # non promoter features
  notPromRFs <- NULL
  for(epi in names(epiPeaks)){
    notPromRFs <- rbind(notPromRFs, 
                        cbind(bt.intersect(
                          a = rfs[rfs[,'type'] != 'promoter',
                                  c('chr', 'core_start', 'core_end', 'type')],
                          b = epiPeaks[[epi]], u = TRUE, wa = TRUE,
                          f = featureOverlap, F = featureOverlap,
                          e = eVal), V5=epi))
  }
  notPromRFs %>%
    mutate(V6 = paste(V4, V1, paste(V2, V3, sep = "-"), 
                      sep = ':')) -> notPromRFs

  # For promoters we request a 1-bp overlapping. Because we extended UPs that
  # overlapped a TSSs window, it could happen that we identified a promoter
  # but then they do not overlap with union peaks 
  rfs %>%
    filter(type == 'promoter') %>%
    dplyr::select(chr, core_start, core_end, type ) -> promsToOverlap
  
  annotatedTSSs %>%
    bt.intersect(promsToOverlap,
                 u = T, wa = T) -> annotatedTSSsInP

  promsInepis <- NULL
  for(epi in names(epiPeaks)){
    promsInEpi <- NULL
    ocp <- epiPeaks[[epi]]
    OCPInTSS <- bt.intersect(ocp, annotatedTSSsInP, wa = T, wb = T)
    if(nrow(OCPInTSS) > 0){
      OCPInTSS %>% 
        group_by(V1,V2,V3) %>% 
        summarise(chr=unique(V1), 
                  start=min(min(V5), V2), 
                  end=max(max(V6), V3), .groups = 'drop') %>% 
        as.data.frame() -> OCPInTSSExtended
      ov <- bt.intersect(a = promsToOverlap, 
                         b = OCPInTSSExtended[,c('chr', 'start', 'end')],
                         u = TRUE, wa = TRUE)
      if(nrow(ov) > 0){
        promsInEpi <- rbind(promsInEpi, cbind(ov, V5=epi)) 
      }
    }
    OCPNotInTSS <- bt.intersect(a = ocp, b = annotatedTSSsInP, wa = T, v=T)
    OCPInP <- bt.intersect(a = promsToOverlap, 
                           b = OCPNotInTSS,
                           u = TRUE, wa = TRUE)
    if(nrow(OCPInP) > 0){
      promsInEpi <- rbind(promsInEpi, 
                         cbind(OCPInP, V5=epi)) 
    }
    promsInepis <- rbind(promsInepis, unique(promsInEpi))
  } 
  promsInepis %>%
    mutate(V6 = paste(V4, V1, paste(V2, V3, sep = "-"),  sep = ':')) %>%
    bind_rows(notPromRFs) %>%
    dplyr::select(-V4) %>%
    rename(chr=V1, start=V2, end=V3, epi=V5, loc=V6)-> rfsInEpiLong
    return(rfsInEpiLong)
}
exportEpigenomeActivity <- function(outPath, species, assembly, rfsGFF, 
                                    epigenomes, rfsInEpiLong){
  date <- paste0(strsplit(as.character(Sys.Date()), 
                          split = "-")[[1]], collapse = "")
  rfsGFF %>%
    filter(type != 'EMAR') -> rfsGFF
  rfsGFF %>%
    select(attributes) %>%
    mutate(stable_ID = str_split(str_split(attributes, 
                                           ";", simplify = TRUE)[,1],
                                 "=", simplify = TRUE)[,2]) %>%
    pull(stable_ID) -> stable_ID
  for(epi_i in epigenomes){
    rfsInEpiLong %>%
      filter(epi == epi_i) %>%
      pull(loc) -> rfInEpi
    
    actLev <- data.frame(stable_ID = stable_ID, status = 'NO EVIDENCE',
                         row.names = rownames(rfsGFF))
    actLev[rfInEpi, 'status'] <- 'EVIDENCE'
    write.csv(actLev, file = paste(outPath,
                                   paste(species, assembly, epi_i,
                                         "Regulatory_activity",
                                         date, "csv", sep = "."),
                                   sep = "/"),
              row.names = F, quote=F)
  }
}

main <- function(parser){
  args <- parser$parse_args()
  species <- args$species
  assembly <- args$assembly
  ESpePre <- args$ESpePre
  startID <- as.numeric(args$startID)
  # read TSSs annotation file
  annotatedTSSsAndInfo <- tryCatch(read.delim(args$annotatedTSSs, header=F, 
                                              sep= "\t"), 
                                   warning = function(w) NULL) 
  stopifnot("File with annotated TSSs doesn't exist" = 
              !is.null(annotatedTSSsAndInfo))
  # extract chromosome names 
  annotatedTSSsAndInfo %>%
    pull(V1) %>% unique() -> chrs
  # read CDSs annotation file
  annotatedCDS <- tryCatch(read.delim(args$annotatedCDS, header =F),
                             warning = function(w) NULL)
  stopifnot("File with annotated CDSs doesn't exist" = 
              !is.null(annotatedCDS))
  path2OCFiles <- args$path2OCFiles
  path2K27acFiles <- args$path2K27acFiles
  path2K4me1Files <- args$path2K4me1Files
  path2K4me3Files <- args$path2K4me3Files
  # if neither k4me1 or k27ac files were specified, then enhancers cannot be 
  # annotated
  if_else(is.null(path2K4me1Files) & is.null(path2K27acFiles), 
          FALSE, TRUE) -> annotateEnhancers
  outPath <- args$outPath
  if(!file.exists(outPath)){
    system(paste('mkdir ', outPath, sep= ''))
  }
  oldAnnotation <- args$oldAnnotation
  if(!is.null(oldAnnotation)){
    # if exist, load previous regulatory features
    previousRegBuild <- tryCatch(read.delim(oldAnnotation, header =F, skip = 1),
                                 warning = function(w) NULL)
    stopifnot("File with previous regulatory annotation doesn't exist" = 
                !is.null(oldAnnotation))
    #Extracting extended region coordinates, feature type and ID
    previousRegBuild %>%
      dplyr::select(V1, V3, V4, V5, V9) %>%
      rename('chr'= V1, 'start' = V4, 'end' = V5, 'type' = V3,
             'attributes' = V9) %>%
      mutate(extended_start = if_else(grepl('start', attributes),
                                      as.numeric(str_split(str_split(attributes, 
                                                                     '_start=', 
                                                                     simplify = TRUE)[,2],
                                                           ';', 
                                                           simplify = T)[,1]),
                                      start)) %>%
      mutate(extended_end = if_else(grepl('end', attributes),
                                    as.numeric(str_split(str_split(attributes,
                                                                   '_end=', 
                                                                   simplify = TRUE)[,2],
                                                         ';',
                                                         simplify = T)[,1]),
                                    end)) %>%
      mutate(ID = str_split(str_split(attributes, 'ID=', 
                                      simplify = TRUE)[,2],
                            ';', simplify = T)[,1]) %>% 
      dplyr::select(chr, extended_start, extended_end, type, ID) -> previousRegBuild
    #Extracting latest ID
    previousRegBuild %>% 
      dplyr::select(ID) %>%
      mutate(ID = str_split(ID, paste(ESpePre, 'R', sep= ''), 
                            simplify = T)[,2]) %>% 
      pull(ID) %>% as.numeric()%>% max() -> oldEndID
    # if there is an old annotation, the new feature IDs will start right after 
    # the last one in the current set
    if(startID <= oldEndID){
      startID <- oldEndID + 1
      print('The startID is set based on old regulatory annotation')
    }
  }else{
    previousRegBuild <- NULL
  }
  ocPattern <- args$ocPattern
  k4me1Pattern <- args$k4me1Pattern
  k27acPattern <- args$k27acPattern
  k4me3Pattern <- args$k4me3Pattern
  TSSwindowUp <- as.numeric(args$TSSwindowUp)
  TSSwindowDown <- as.numeric(args$TSSwindowDown)
  coreWUp <- as.numeric(args$coreWUp)
  coreWDown <- as.numeric(args$coreWDown)
  featureOverlap <- as.numeric(args$featureOverlap)
  eVal <- FALSE
  if(!is.null(featureOverlap)){
    featureOverlap <- as.numeric(featureOverlap)
    eVal <- TRUE 
  }
  histoneOverlap <- as.numeric(args$histoneOverlap)
  OCRsHistoneOverlap <- as.numeric(args$OCRsHistoneOverlap)
  cOCRsCDSOverlap <- as.numeric(args$cOCRsCDSOverlap)
  maxDistOCR <- as.numeric(args$maxDistOCR)
  maxDistH <- as.numeric(args$maxDistH)
  maxPromExtLen <- args$maxPromExtLen
  if(maxPromExtLen != ''){
    maxPromExtLen <- as.numeric(maxPromExtLen)
  }else{
    maxPromExtLen <- NULL
  }
  minOCRLen <- args$minOCRLen
  if(minOCRLen!= ''){
    minOCRLen <- as.numeric(minOCRLen)
  }else{
    minOCRLen <- NULL
  }
  maxEnhancerLen <- args$maxEnhancerLen
  if(maxEnhancerLen!= ''){
    maxEnhancerLen <- as.numeric(maxEnhancerLen)
  }else{
    maxEnhancerLen <- NULL
  }
  maxOCRLen <- args$maxOCRLen
  if(!is.null(maxOCRLen)){
    maxOCRLen<- as.numeric(maxOCRLen)
  }
  exportActivity <- as.logical(args$exportActivity)
  path2EpigenomesOCFiles <- args$path2EpigenomesOCFiles
  if(exportActivity & is.null(path2EpigenomesOCFiles)){
    path2EpigenomesOCFiles <- path2OCFiles
  }
  epiPattern <- args$epiPattern
  
  # The regulatory features will be defined by the following steps
  
  # 1. Load all ATAC-seq/DNase-seq, H3K4me1 and H3K27ac peaks identified across
  # all epigenomes, keeping just peaks in chromosomes of interest. 
  
  # OC peaks
  ocBedFiles <- list.files(path = path2OCFiles, pattern = ocPattern)
  OCPeaks <- loadPeaks(path2OCFiles, ocBedFiles, chrs)
  print(paste("Number of loaded open chromatin bed files: ", 
              length(ocBedFiles), sep = ""))
  # K4me1 peaks
  if(!is.null(path2K4me1Files)){
    k4me1BedFiles <- list.files(path = path2K4me1Files, 
                                pattern = k4me1Pattern)
    h3k4me1Peaks <- loadPeaks(path2K4me1Files, k4me1BedFiles, chrs)
    print(paste("Number of loaded H3K4me1 bed files: ", length(k4me1BedFiles), 
                sep = ""))
  }else{
    h3k4me1Peaks <- NULL
  }
  # K27ac peaks
  if(!is.null(path2K27acFiles)){
    k27acBedFiles <- list.files(path = path2K27acFiles, 
                                pattern = k27acPattern)
    h3k27acPeaks <- loadPeaks(path2K27acFiles, k27acBedFiles, chrs)
    print(paste("Number of loaded H3K27ac bed files: ", length(k27acBedFiles), 
                sep = ""))
  }else{
    h3k27acPeaks <- NULL
  }
  # K4me3 peaks
  if(!is.null(path2K4me3Files)){
    k4me3BedFiles <- list.files(path = path2K4me3Files, 
                                pattern = k4me3Pattern)
    h3k4me3Peaks <- loadPeaks(path2K4me3Files, k4me3BedFiles, chrs)
    print(paste("Number of loaded H3K4me3 bed files: ", length(k4me3BedFiles), 
                sep = ""))
  }else{
    h3k4me3Peaks <- NULL
  }
  
  # 2. Combine all ATAC-seq/DNAse-seq peaks and merge the overlapping ones to 
  #    obtain a set of unique peaks (UPs).
  UPs <- defineUnionPeaks(OCPeaks)
  print(paste("Number of open chromatin union peaks: ", nrow(UPs), sep = ""))
  
  # 3. Identify UPs overlapping with e! TSSs windows for splitting the set of UPs 
  # into UPs-TSSs and UPs-noTSSs.
  # TSSs file has chr,start,end,strand. Window start and end will be added
  annotatedTSSs <- addingTSSsInfo(annotatedTSSsAndInfo, TSSwindowUp,
                                  TSSwindowDown, coreWUp, coreWDown)
  
  UPsTSS <- definingUPsInTSS(UPs, annotatedTSSs)
  print(paste("Number of open chromatin union peaks overlapping TSSs: ",
              nrow(UPsTSS), sep = ""))
  
  # 4. Defining UPs-TSSs as promoters considering the updated start and end and
  # the possibility of UPs overlapping more than 1 TSSs. In this case, if core
  # regions of potential promoters overlap, then the promoter is unique and its
  # boundaries are defined by combining the updated boundaries. If not, then 
  # the UPs-TSS is split having one promoter per TSSs. All UPs-TSSs region not
  # covered by updated boundaries are defined as candidate OCRs
  
  rfs <- promotersIdentification(UPsTSS, TSSwindowDown, maxPromExtLen)
  promoters <- rfs[['promoters']]
  print(paste("Number of identified promoters: ", nrow(promoters), sep = ""))  
  
  # 5. Consider cOCRs and UPs-noTSSs (merged nearby peaks up to maxDistOCR apart)
  #   as merged candidate Open Chromatin Regions (mcOCRs)
  
  cOCRs <- rfs[['cOCRs']] # cOCRs do not overlapping UPs in TSS
  upsIncOCRs <- rfs[['upsIncOCRs']] # UPs in TSS containing cOCRs
  if(nrow(upsIncOCRs) > 0 ){
    upsIncOCRs %>% mutate(V1 = as.character(V1)) -> upsIncOCRs
  }
  # remove the regions of the UPs that overlapp TSSs
  UPs %>% 
    bt.intersect(annotatedTSSs, v = T, wa = T) -> UPsnoTSS
  UPsnoTSS %>% mutate(V1 = as.character(V1)) -> UPsnoTSS
  UPsnoTSS %>%
    bind_rows(upsIncOCRs) -> UPsnoTSS
  # UPsnoTSS <- rbind(UPsnoTSS, upsIncOCRs)
  names(UPsnoTSS) <- names(cOCRs)
  # UPsnoTSS %>% mutate(chr = as.character(chr)) -> UPsnoTSS
  mcOCRs <- OCRsIdentification(UPsnoTSS, cOCRs, promoters, maxDistOCR,
                               minOCRLen)
  # 6. Define merged union peaks for open chromatin peaks and histone marks
  if(!is.null(h3k4me1Peaks) | !is.null(h3k27acPeaks)){
    ocUPs <- defineUnionPeaks(OCPeaks, d = maxDistOCR)
    ocUPs %>%
      bt.subtract(promoters[, c('chr', "extended_start", 
                                "extended_end")]) -> ocUPs
  }
  if(!is.null(h3k4me1Peaks)){
    H3K4me1UPs <- defineUnionPeaks(h3k4me1Peaks, d = maxDistH)
    print(paste("Number of merged union H3K4me1 peaks: ", 
                nrow(H3K4me1UPs), sep = ""))
  }
  if(!is.null(h3k27acPeaks)){
    H3K27acUPs <- defineUnionPeaks(h3k27acPeaks, d = maxDistH)
    print(paste("Number of merged union H3K27acUPs peaks: ", 
                nrow(H3K27acUPs), sep = ""))  
  }
  # 7. Identify H3K4me1 UPs and H3K27ac UPs overlapping with OC union peaks.
  if(!is.null(h3k4me1Peaks)){
    H3K4me1OCUPs <- bt.intersect(a = H3K4me1UPs, b = ocUPs, e = T, 
                                 f = histoneOverlap, F = histoneOverlap, 
                                 wa = T, u = T)
    names(H3K4me1OCUPs) <- c('chr', 'start', 'end')
    print(paste("Number of H3K4me1-OC union peaks: ", nrow(H3K4me1OCUPs),
                sep = ""))
  }else{
    H3K4me1OCUPs <- NULL
  }
  if(!is.null(h3k27acPeaks)){
    H3K27acOCUPs <- bt.intersect(a = H3K27acUPs, b = ocUPs, e = T, 
                                 f = histoneOverlap, F = histoneOverlap, 
                                 wa = T, u = T)
    names(H3K27acOCUPs) <- c('chr', 'start', 'end')
    print(paste("Number of H3K27ac-OC union peaks: ", nrow(H3K27acOCUPs),
                sep = ""))
  }else{
    H3K27acOCUPs <- NULL
  }
  # 8. If histone peaks are available, identify enhancers by overlapping
  # mcOCRs with H3K4me1OCUPs and H3K27acOCUPs peaks 
  if(annotateEnhancers){
    rfs <- enhancersIdentification(mcOCRs, minOCRLen, H3K4me1OCUPs, 
                                   H3K27acOCUPs, annotatedCDS, cOCRsCDSOverlap,
                                   OCRsHistoneOverlap, maxEnhancerLen)
    enhancers <- rfs[['enhancers']]
    rownames(enhancers) <- paste(enhancers[,'chr'], 
                                 paste(enhancers[,'extended_start'], 
                                       enhancers[,'extended_end'], sep = '-'),
                                 sep = ':')
    print(paste("Number of identified enhancers: ", nrow(enhancers),
              sep = ""))
  }
  # 9. Final set of OCRs
  if(annotateEnhancers){
    mcOCRs <- rfs[['mcOCRs']]
  }
  if(!is.null(maxOCRLen)){
    mcOCRs %>%
      filter(extended_end - extended_start <= maxOCRLen) -> mcOCRs
  }
  rownames(mcOCRs) <- paste(mcOCRs[,'chr'],
                            paste(mcOCRs[,'extended_start'],
                                  mcOCRs[,'extended_end'], sep = '-'),
                            sep = ':')
  print(paste("Number of identified open chromatin regions: ",
              nrow(mcOCRs), sep = ""))  
  # 10. EMARs identification
  if(!is.null(h3k4me1Peaks) | !is.null(h3k27acPeaks) | !is.null(h3k4me3Peaks)){
    ocUPs <- defineUnionPeaks(OCPeaks, d = maxDistOCR)
    #for considering promoters and the other features independently
    ocUPs %>% 
      bt.subtract(promoters[,c('chr', 'extended_start',
                               'extended_end')]) %>%
      bind_rows(bt.intersect(ocUPs, promoters[,c('chr', 'extended_start',
                                                 'extended_end')])) -> ocUPs
    if(!is.null(h3k4me1Peaks)){
      H3K4me1OCUPs <- bt.intersect(a = H3K4me1UPs, b = ocUPs, e = T, 
                                   f = histoneOverlap, F = histoneOverlap, 
                                   wa = T, u = T)
      names(H3K4me1OCUPs) <- c('chr', 'start', 'end')
      print(paste("Number of H3K4me1-OC union peaks: ", nrow(H3K4me1OCUPs),
                  sep = ""))
    }else{
      H3K4me1OCUPs <- NULL
    }
    if(!is.null(h3k27acPeaks)){
      H3K27acOCUPs <- bt.intersect(a = H3K27acUPs, b = ocUPs, e = T, 
                                   f = histoneOverlap, F = histoneOverlap, 
                                   wa = T, u = T)
      names(H3K27acOCUPs) <- c('chr', 'start', 'end')
      print(paste("Number of H3K27ac-OC union peaks: ", nrow(H3K27acOCUPs),
                  sep = ""))
    }else{
      H3K27acOCUPs <- NULL
    }
    if(!is.null(h3k4me3Peaks)){
      H3K4me3UPs <- defineUnionPeaks(h3k4me3Peaks, d = maxDistH)
      print(paste("Number of merged union H3K4me3UPs peaks: ", 
                  nrow(H3K4me3UPs), sep = ""))  
      H3K4me3OCUPs <- bt.intersect(a = H3K4me3UPs, b = ocUPs, e = T,
                                   f = histoneOverlap, F = histoneOverlap,
                                   wa = T, u = T)
      names(H3K4me3OCUPs) <- c('chr', 'start', 'end')
      print(paste("Number of H3K4me3-OC union peaks: ", nrow(H3K4me3OCUPs),
                  sep = ""))
    }else{
      H3K4me3OCUPs <- NULL
    }
    emars <- emarsIdentification(ocUPs, H3K4me1OCUPs, H3K27acOCUPs,
                                 H3K4me3OCUPs, OCRsHistoneOverlap, minOCRLen)
    print(paste("Number of identified EMARs: ", nrow(emars), sep = ""))  
  }else{
    emars <- NULL
  }
  # 11. Generate RFs gff table
  # RFs will be combined in a single object, adding the regulatory feature type
  promoters %>% 
    mutate(chr = as.character(chr), type = 'promoter') -> rfs
  
  if(annotateEnhancers){
    rfs %>%
      bind_rows(enhancers %>% mutate(chr = as.character(chr),
                                     type = 'enhancer')) -> rfs
  }
  rfs %>%
    bind_rows(mcOCRs %>% 
                mutate(chr = as.character(chr), 
                       type = 'open_chromatin_region')) -> rfs
    
  if(!is.null(emars)){
    rfs %>%
      bind_rows(emars %>% 
                  mutate(chr = as.character(chr), 
                         type = 'EMAR')) -> rfs
  }
  names(rfs) <- c(names(promoters), 'type')
  if(!is.null(previousRegBuild)){
    featuresToMerge <- c('promoter', 'open_chromatin_region')
    if(annotateEnhancers) featuresToMerge <- c(featuresToMerge, 'enhancer')
    if(!is.null(emars)) featuresToMerge <- c(featuresToMerge, 'EMAR')
  }
  RFsGFF <- createRFsGFF(rfs, chrs, ESpePre, startID, previousRegBuild, 
                         annotatedTSSsAndInfo, featuresToMerge = featuresToMerge, 
                         source = 'Ensembl')
  date <- paste0(strsplit(as.character(Sys.Date()), 
                          split = "-")[[1]], collapse = "")
  gff3h <- "##gff-version 3"
  RFFile <- paste(outPath, paste(species, assembly, "regulatory_features", 
                                 date, "gff3", sep = "."), sep = "/")
  RFFileCon <- file(RFFile)
  writeLines(text = gff3h, con = RFFileCon)
  close(RFFileCon)
  write.table(RFsGFF, file = RFFile, 
              col.names = F, row.names = F, sep = "\t", quote = F,
              append = T)
  if(exportActivity){
    epiBedFiles <- list.files(path = path2EpigenomesOCFiles, 
                             pattern = epiPattern)
    epiPeaks <- loadEpiPeaks(path2OCFiles, epiBedFiles, chrs, assembly)
    print(paste("Number of epigenomes: ", length(epiPeaks), sep = ""))
    rfsInEpiLong <- rfsInEpigenomes(rfs, epiPeaks, annotatedTSSs, 
                                    featureOverlap, eVal)
    exportEpigenomeActivity(outPath, species, assembly, RFsGFF, names(epiPeaks),
                            rfsInEpiLong)  
  }
  
}

parser <- parsingArgs()
main(parser)
