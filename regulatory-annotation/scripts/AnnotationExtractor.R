# Script: AnnotationExtractor.R
# Author: Gabriela Merino - Ensembl Regulation
# Date: 10-2025
# Description: Extract TSSs and CDSs annotation from Ensembl/GENCODE GFF3 files

suppressPackageStartupMessages({
  library(dplyr)
  library(argparse)
  library(stringr)
  library(bedtoolsr)
})

options(scipen = 999)

parsingArgs <- function(){
  parser <- ArgumentParser(prog = 'Rscript AnnotationExtractor.R',
                           description = 'Extract TSSs and CDSs annotation from Ensembl GFF3 file')
  parser$add_argument('chrs', metavar = 'ChrFile', type = "character",
                      nargs = 1, help = 'Full path to the file listing the 
                      canonical chromosomes / chromosomes of interest (one per line).')
  parser$add_argument('assemblyAnnotation', metavar = 'annotationGFF3',
                      type = "character", nargs = 1,
                      help = 'Full path to the Ensembl gff3 file containing 
                      gene annotation.')
  parser$add_argument('isGENCODE', metavar = 'isGENCODE', type = "logical", 
                      nargs = 1, help = 'Logical indicating if the genome 
                      annotation file has been produced by GENCODE.')
  parser$add_argument('--transcriptsBiotypes', dest = 'transcriptsBiotypes', 
                      action = 'store', default = NULL, 
                      help = 'Comma-separated list of transcript biotypes to 
                      keep. Default value is `NULL` indicating that all biotypes
                      should be kept.')
  parser$add_argument('--flagsToKeep', dest = 'flagsToKeep', 
                      action = 'store', default = NULL, 
                      help = 'Comma-separated list of transcript flags 
                      to be keep. Default value is `NULL` indicating that all 
                      transcripts with any flag are kept.')
  parser$add_argument('--flagsToRemove', dest = 'flagsToRemove', 
                      action = 'store', default = NULL, 
                      help = 'Comma-separated list of flags 
                      indicating transcripts having them should be removed. 
                      Default value is NULL, all transcripts are kept.')
  return(parser)
}

extractTSSs <- function(genomeAnnot, transcriptsBiotypes, isGENCODE, flagsToKeep,
                        flagsToRemove){
  # keep only transcript-level annotation
  genomeAnnot %>%
    filter(type == 'transcript') %>%
    select(-type) -> transcriptsAnnot
  # extract transcript ID and biotype
  if(isGENCODE){
    transcriptsAnnot %>%
      mutate(ID = str_split(str_split(attributes, 'ID=', 
                                      simplify = T)[,2],
                            '[.]', simplify = T)[,1],
             biotype = str_split(str_split(attributes, 'transcript_type=', 
                                           simplify = T)[,2],
                                 ';', simplify = T)[,1]) -> transcriptsAnnot
  }else{
    transcriptsAnnot %>%
      mutate(ID = str_split(str_split(attributes, 'ID=transcript:', 
                                      simplify = T)[,2],
                            ';', simplify = T)[,1],
             biotype = str_split(str_split(attributes, 'biotype=', 
                                           simplify = T)[,2],
                                 ';', simplify = T)[,1]) -> transcriptsAnnot
  }
  # Filter transcripts based on the specified biotypes
  if(!is.null(transcriptsBiotypes)){
    transcriptsAnnot %>%
      filter(biotype %in% transcriptsBiotypes) -> transcriptsAnnot
  }
  # Filter transcripts based on the specified flags
  if(!is.null(flagsToKeep)){
    transcriptsAnnot %>%
      filter(grepl(flagsToKeep, attributes)) -> transcriptsAnnot
  }
  if(!is.null(flagsToRemove)){
    transcriptsAnnot %>%
      filter(!grepl(flagsToRemove, attributes)) -> transcriptsAnnot
  }
  # Defining TSSs coordinates
  transcriptsAnnot %>%
    mutate(end = if_else(strand == '+', start + 1, end)) %>%
    mutate(start = if_else(strand == '-', end -1, start)) -> transcriptsAnnot
  # Assigning gene ID to each TSS
  if(isGENCODE){
    transcriptsAnnot %>%
      mutate(gene = str_split(str_split(attributes, 'Parent=', 
                                        simplify = T)[,2],
                              '[.]', simplify = T)[,1]) -> transcriptsAnnot
  }else{
    transcriptsAnnot %>%
      mutate(gene = str_split(str_split(attributes, 'Parent=gene:', 
                                        simplify = T)[,2],
                              ';', simplify = T)[,1]) -> transcriptsAnnot
  }
  transcriptsAnnot %>%
    select(-attributes) -> transcriptsAnnot
  return(transcriptsAnnot)
}
addGeneInfo <- function(genomeAnnot, transcriptsAnnot, isGENCODE){
  # keep only gene-level annotation
  genomeAnnot %>%
    filter(type == 'gene') %>%
    select(chr, start, end, strand, attributes) -> geneAnnot
  # extract gene ID, biotype, and name
  if(isGENCODE){
    geneAnnot %>%
      mutate(ID = str_split(str_split(attributes, 'ID=', simplify = T)[,2],
                            '[.]', simplify = T)[,1],
             gene_biotype = str_split(str_split(attributes, 'gene_type=', 
                                                simplify = T)[,2],
                                      ';', simplify = T)[,1],
             gene_name = str_split(str_split(attributes, 'gene_name=', 
                                                simplify = T)[,2],
                                      ';', simplify = T)[,1]) -> geneAnnot
  }else{
    geneAnnot %>%
      mutate(ID = str_split(str_split(attributes, 'ID=gene:', simplify = T)[,2],
                            ';', simplify = T)[,1],
             gene_biotype = str_split(str_split(attributes, 'biotype=', 
                                                simplify = T)[,2],
                                      ';', simplify = T)[,1]) -> geneAnnot
  }
  # Add gene metadata to transcripts table
  transcriptsAnnot %>%
    left_join(geneAnnot %>% 
                select(-c(chr, start, end, strand, attributes)) %>% 
                rename(gene = ID)) %>%
    select(-biotype) -> transcriptsAnnot
  return(transcriptsAnnot)
}
  
extractCDSs <- function(genomeAnnot, transcriptsAnnot, isGENCODE){
  # keep only CDS-level annotation
  genomeAnnot %>%
    filter(type == 'CDS') %>%
    select(chr, start, end, attributes) -> CDSAnnot
  # extract transcript ID
  if(isGENCODE){
    CDSAnnot %>%
      mutate(transcript = str_split(str_split(attributes, 'Parent=',
                                              simplify = T)[,2], '[.]',
                                    simplify = T)[,1]) -> CDSAnnot
  }else{
    CDSAnnot %>%
      mutate(transcript = str_split(str_split(attributes, 'Parent=transcript:',
                                              simplify = T)[,2],
                                    ';', simplify = T)[,1]) -> CDSAnnot
  }
  # keep CDS from transcripts of interest and merge overlapping regions
  CDSAnnot %>%
    select(-attributes) %>%
    filter(transcript %in% (transcriptsAnnot %>% pull(ID))) %>%
    select(chr, start, end) %>%
    bt.sort() %>%
    bt.merge() -> CDSAnnot
  return(CDSAnnot)
}

main <- function(parser){
  args <- parser$parse_args()
  # Read chromosomes name file
  chrs <- tryCatch(read.delim(args$chrs, header = F),
                   warning = function(w) NULL)[,1]
  stopifnot("File with chromosomes names doesn't exist" = 
              !is.null(chrs))
  assemblyAnnotation <- args$assemblyAnnotation
  # Read annotation file
  genomeAnnot <- tryCatch(read.delim(assemblyAnnotation, header = F, 
                                     sep= "\t", comment.char = '#'),
                          warning = function(w) NULL) 
  names(genomeAnnot) <- c('chr', 'source', 'type', 'start', 'end', 'score',
                          'strand', 'phase', 'attributes')
  stopifnot("Annotation file doesn't exist" = !is.null(genomeAnnot))
  isGENCODE <- as.logical(args$isGENCODE)
  # For compatibility with fasta format
  if(isGENCODE){
    genomeAnnot %>%
      mutate(chr = str_split(chr, 'chr', simplify = T)[,2]) -> genomeAnnot 
  }
  # Filter to keep chromosomes of interest
  genomeAnnot %>%
    filter(chr %in% chrs) %>%
    mutate(chr = factor(chr, levels = chrs)) %>%
    arrange(chr, start, end) -> genomeAnnot
  
  # Extract transcripts biotype of interest
  if(!is.null(args$transcriptsBiotypes)){
    transcriptsBiotypes <- strsplit(args$transcriptsBiotypes, split = ',')[[1]]
  }
  # Extract transcripts flags to keep
  if(!is.null(args$flagsToKeep)){
    flagsToKeep <- strsplit(args$flagsToKeep, split = ',')[[1]]
  }else{
    flagsToKeep <- NULL
  }
  # Extract transcripts flags to remove
  if(!is.null(args$flagsToRemove)){
    flagsToRemove <- strsplit(args$flagsToRemove, split = ',')[[1]]
  }else{
    flagsToRemove <- NULL
  }
  # Selecting columns of interest and modifying coordinates to be in zero-based
  # (bed format)
  genomeAnnot %>%
    select(chr, start, end, strand, type, attributes) %>%
    mutate(start = start - 1) -> genomeAnnot 
  # For Ensembl annotation, extract feature type from attributes.  
  if(!isGENCODE){  
    genomeAnnot %>%
      mutate(type = str_split(str_split(attributes, 'ID=', simplify = T)[,2],
                              ':', simplify = T)[,1]) -> genomeAnnot
  }
  # Keep gene, transcript and CDS features
  genomeAnnot  %>%
    filter(type %in% c('gene', 'transcript', 'CDS')) -> genomeAnnot
  # Extract TSSs annotation
  transcriptsAnnot <- extractTSSs(genomeAnnot, transcriptsBiotypes, isGENCODE,
                                  flagsToKeep, flagsToRemove)
  # Add gene information 
  transcriptsAnnot <- addGeneInfo(genomeAnnot, transcriptsAnnot, isGENCODE)
  CDSAnnot <- extractCDSs(genomeAnnot, transcriptsAnnot, isGENCODE)
  # Export TSS and CDS annotation files
  prefix <- strsplit(assemblyAnnotation, '.gff3')[[1]][1]
  write.table(transcriptsAnnot %>% select(-ID) %>% unique(),
              paste(prefix, '_TSS.bed', sep =''), col.names = FALSE, 
              row.names = FALSE, quote = FALSE, sep = '\t')
  
  write.table(CDSAnnot, paste(prefix, '_CDS.bed', sep =''),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}

parser <- parsingArgs()
main(parser)
