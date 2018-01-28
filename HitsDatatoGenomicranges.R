AnnotatingFromHitlist <- function(Hitlist, PatternList){
  require(GenomicRanges); require(tidyverse)
  Hits <- read.csv(Hitlist, header=T, sep ="\t")
  #create Granges Object with DmelCoordinates
  a <- separate(Hits, genome_pos_mel, into = c("seqnames", "Start", "Stop"))
  a.df <- as.data.frame(seqnames = seqnames, sTarT = Start, sToP = Stop)
  GrangesObject <<- makeGRangesFromDataFrame(a.df)
  ## Add Coordinates For other species as MetadataColumn
  mcols(GrangesObject)$yakCoord <- Hits$genome_pos_yak
  mcols(GrangesObject)$pseCoord <- Hits$genome_pos_pse
  mcols(GrangesObject)$virCoord <- Hits$genome_pos_vir
  GrangesObject
  ##PatternList 
  Pattern <- read.csv(PatternList, header=T, sep ="\t")
  PatternA <- Pattern[,c("VTID","Chrosome","Start","End", "stg4_6")]
  #taking only stg4_6 only gets 
  PatternA <- na.omit(PatternA)
  ##obtaining the most expressed pattern of each contruction the ones with equal exp are mainteined
  df_lines <- PatternA %>% 
    separate_rows(stg4_6, sep = '\\|') %>% 
    separate(pattern, c('Pattern', 'LevelofExp'), sep = ';', convert = TRUE)
  clean.list.df   <- df_lines %>% group_by(VTID) %>% top_n(1, value)
  colnames(clean.list.df)[2] <- "seqnames"
  colnames(clean.list.df)[4] <- "sToP"
  clean.list.gr <- GRanges(clean.list.df)
  
}