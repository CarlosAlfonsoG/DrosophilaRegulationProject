makePhasConsBedGranges <- function(fileA, fileB, TF, Tissue){
  PhasConsValuesfromBed <- as.data.frame(read.table(fileA, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")) 
  PhasConsValuesfromBed.df <- data.frame(seqnames=PhasConsValuesfromBed$V1, sTarT = PhasConsValuesfromBed$V2, sToP = PhasConsValuesfromBed$V3)
  PhasConsValuesfromBed.gr <- makeGRangesFromDataFrame(PhasConsValuesfromBed.df)
  score(PhasConsValuesfromBed.gr) <- PhasConsValuesfromBed$V5
  PeakAnno <- read.csv(fileB, sep=",", header = T)
  PeakAnno$Pattern <- paste(Tissue, PeakAnno$Pattern)
  PeakAnno$TF <- paste(Tissue, PeakAnno$TF)
  PeakAnno.df <- data.frame(seqnames=PeakAnno$seqnames, sTarT=PeakAnno$start, sTop=PeakAnno$end)
  PeakAnno.gr <- makeGRangesFromDataFrame(PeakAnno.df)
  mcols(PeakAnno.gr)$Pattern <- PeakAnno$Pattern 
  mcols(PeakAnno.gr)$annotaton <- PeakAnno$annotation
  mcols(PeakAnno.gr)$distanceToTSS <- PeakAnno$distanceToTSS
  mcols(PeakAnno.gr)$geneId <- PeakAnno$geneId
  mcols(PeakAnno.gr)$transcriptId <- PeakAnno$transcriptId
  mcols(PeakAnno.gr)$TF <- PeakAnno$TF
  ranges <- subsetByOverlaps(PeakAnno.gr, PhasConsValuesfromBed.gr)
  hits <-findOverlaps(PeakAnno.gr, PhasConsValuesfromBed.gr)
  idx <- unique(subjectHits(hits))  
  values <- DataFrame(PhasConsScore=PhasConsValuesfromBed.gr$score[idx])
  mcols(ranges) <- c(mcols(ranges), values)
  varname <- paste(TF, Tissue, sep = "_")
  Annotation <<- as.data.frame(ranges)
  write.csv(Annotation, varname)
  
}
