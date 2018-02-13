AnnotateGrangesWithPatternMetada <- function(Hitlist, PatternList, GeneName){
  options(verbose = F)
  require(GenomicRanges); require(tidyverse);require(GenomicRanges); require(gsubfn):require(stringi)
  Hits <- read.table(Hitlist, header=T, sep ="\t")
  #create Granges Object with DmelCoordinates
  a <- separate(Hits, genome_pos_mel, into = c("seqnames", "Start", "Stop"))
  a.df <- data.frame(seqnames = a$seqnames, sTarT = a$Start, sToP = a$Stop)
  GrangesObject <- makeGRangesFromDataFrame(a.df)
  ## Add Coordinates For other species as MetadataColumn
  mcols(GrangesObject)$yakCoord <- Hits$genome_pos_yak
  mcols(GrangesObject)$pseCoord <- Hits$genome_pos_pse
  mcols(GrangesObject)$virCoord <- Hits$genome_pos_vir
  GrangesObject
  ##PatternList 
  Pattern <- read.csv(PatternList, header=T, sep =",")
  PatternA <- Pattern[,c("VTID","Chrosome","Start","End", "stg4_6")]
  #taking only stg4_6 only gets 
  PatternA <- na.omit(PatternA)
  ##obtaining the most expressed pattern of each contruction the ones with equal exp are mainteined
  df_lines <- PatternA %>% 
    separate_rows(stg4_6, sep = '\\|') %>% 
    separate(stg4_6, c('Pattern', 'LevelofExp'), sep = ';', convert = TRUE)
  clean.list.df   <- df_lines %>% group_by(VTID) %>% top_n(1, LevelofExp)
  colnames(clean.list.df)[2] <- "seqnames"
  colnames(clean.list.df)[4] <- "sToP"
  GRanges.Patterns <- GRanges(clean.list.df)
  GRanges.Patterns
  Gr1 <- GrangesObject
  Gr2 <- GRanges.Patterns
  ranges <-subsetByOverlaps(Gr1, Gr2)
  hits <- findOverlaps(Gr1, Gr2)
  Pattern <- CharacterList(split(Gr2$Pattern[subjectHits(hits)], queryHits(hits)))
  ReporterIntensity <- CharacterList(split(Gr2$LevelofExp[subjectHits(hits)], queryHits(hits)))
  VTID <- CharacterList(split(Gr2$VTID[subjectHits(hits)], queryHits(hits)))
  mcols(ranges) <- DataFrame(mcols(ranges), Pattern, ReporterIntensity,VTID)
  AnnotatedGrangesWithPattern <- ranges
  AnnotatedGrangesWithPattern
  #CleanData and generate Duplicate Rows for pattern that match several coordinates 
  AnnotatedGrangesWithPattern.df <- as.data.frame(AnnotatedGrangesWithPattern)
  AnnotatedGrangesWithPattern.df$Pattern <- as.character(AnnotatedGrangesWithPattern.df$Pattern)
  AnnotatedGrangesWithPattern.df$Pattern <- gsub("\\(", "", AnnotatedGrangesWithPattern.df$Pattern)
  AnnotatedGrangesWithPattern.df$Pattern <- gsub("\\)", "", AnnotatedGrangesWithPattern.df$Pattern)
  AnnotatedGrangesWithPattern.df$Pattern <- gsub('"', "", AnnotatedGrangesWithPattern.df$Pattern)
  AnnotatedGrangesWithPattern.df$Pattern <- gsub("^c", "", AnnotatedGrangesWithPattern.df$Pattern)
  AnnotatedGrangesWithPattern.df$ID <- seq.int(nrow(AnnotatedGrangesWithPattern.df))
  df <- AnnotatedGrangesWithPattern.df
  s <- strsplit(df$Pattern, split = ",")
  e <- data.frame(ID = rep(df$ID, sapply(s, length)), Pattern = unlist(s))
  colnames(e)[1] <- "ID"
  f <- merge(df, e, by="ID", all.y = T)
  f$ReporterIntensity <- as.character(f$ReporterIntensity)
  f$ReporterIntensity <- stri_extract_first_regex(f$ReporterIntensity, "[0-9]")
  DataAnnnobyPattern <- f
  DataAnnnobyPattern$VTID <- NULL
  DataAnnnobyPattern$Pattern.x <- NULL
  ## Create Granges Objects and subset by specific patterns 
  TF.GRangesAndPattern  <<- GRanges(DataAnnnobyPattern)
  Ectoderm <- DataAnnnobyPattern[grep("ecto", DataAnnnobyPattern$Pattern.y), ]
  Endoderm <- DataAnnnobyPattern[grep("endo", DataAnnnobyPattern$Pattern.y), ]
  Mesoderm <- DataAnnnobyPattern[grep("meso", DataAnnnobyPattern$Pattern.y), ]
  Ubiquitous <- DataAnnnobyPattern[grep("ubiquitous", DataAnnnobyPattern$Pattern.y), ]
  GAP <- DataAnnnobyPattern[grep("gap", DataAnnnobyPattern$Pattern.y), ]
  AP <- DataAnnnobyPattern[grep("A-P", DataAnnnobyPattern$Pattern.y), ]
  ## create bed files for each and directories 
  newhome <- c("Results")
  dir.create(newhome)
  setwd(newhome)
  write.table(Ectoderm, file=paste0(GeneName, "_Ectoderm.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(Endoderm, file=paste0(GeneName, "_Endoderm.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(Mesoderm, file=paste0(GeneName, "_MesodermHits.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(Ubiquitous, file=paste0(GeneName, "_UbiquitousHits.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(AP, file=paste0(GeneName, "_A-P.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(GAP, file=paste0(GeneName, "_gap.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}
