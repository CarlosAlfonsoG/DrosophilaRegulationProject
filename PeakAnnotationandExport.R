##Anotating sorted peaksssssss ### 
library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(BSgenome.Dmelanogaster.UCSC.dm3)
###FirstCreate Granges Objects ## 
Range.df <- data.frame(seqnames = df$chrnames, StarT = df$Coord1, sToP = df$Coord2) 
GRangeObjt <- makeGRangesFromDataFrame(Range.df)
###Download the genome to annotate the regions of interest using GenomicFeatures library##
dm3 <- makeTxDbFromUCSC(genome = "dm3", tablename ="ensGene")
# To determine available list of a given genome use:  supportedUCSCtables(genome = "genomeofinterest")
## To AnnotatePeaks just use: 
peakAnnoUbiquitousHun <- annotatePeak(HunUbiquitous , tssRegion=c(-3000, 3000), TxDb=dm3)
## Pieplot
plotAnnoPie(peakAnnoMesodermHun)
##Export coordinates to FASTA files ## 
seq = BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm3, regionofinterest)
Biostrings::writeXStringSet(seq, "../FastaFiles/BicoidEctoderm.out.fasta")


