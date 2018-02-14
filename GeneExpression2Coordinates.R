###This function works specifically with EisenMetadata and GRangesObjects To annotate them with nearby genes using ChipSeeker ## 
## Need upgrade for the following things, be more generalized for other usages, GeneExpression value is useless, create a directory and be able to run with several BED files at a time ##
GeneExpression2Coordinates <- function(BED, GeneExpresion, GeneName){
    require(GenomicRanges):require(ChIPseeker):require(GenomicFeatures)
    data <- read.table(BED,header=F)
    colnames(data) <- c('ID','chr','start', 'end','size','strand', 'Dyak', 'Dpse', 'Dvir', 'PatternInt', 'Pattern')
    GRangesChipSeeker <- makeGRangesFromDataFrame(data)
    mcols(GRangesChipSeeker)$Dyak <- data$Dyak
    mcols(GRangesChipSeeker)$Dpse <- data$Dpse
    mcols(GRangesChipSeeker)$Dvir <- data$Dvir
    mcols(GRangesChipSeeker)$PatternExp <- data$PatternInt
    mcols(GRangesChipSeeker)$Pattern <- data$Pattern
    if(  !exists("dm3") ) {
      dm3 <- makeTxDbFromUCSC(genome = "dm3", tablename ="ensGene")
    }
    peakAnnot <- annotatePeak(GRangesChipSeeker, tssRegion=c(-3000, 3000), TxDb=dm3)
    data.frame <<- as.data.frame(peakAnnot)
    colnames(data.frame)[17] <- "Gene"
    data.frame$annotation <- gsub(",","",data.frame$annotation)
    Norm_RNAseq <- read.csv("/Users/carlosalfonso/Projects/EnhancerDivergence/Data/Tables/csv/Norm_Gene_exp_all_species.csv", sep = ",", header = T)
    colnames(Norm_RNAseq)[1] <- "Gene"
    data.frame  <- merge(x = data.frame, y = Norm_RNAseq, by = "Gene", all.x = TRUE)
    Ectoderm <- data.frame[grep("ecto", data.frame$Pattern), ]
    Endoderm <- data.frame[grep("endo", data.frame$Pattern), ]
    Mesoderm <- data.frame[grep("meso", data.frame$Pattern), ]
    Ubiquitous <- data.frame[grep("ubiquitous", data.frame$Pattern), ]
    GAP <- data.frame[grep("gap", data.frame$Pattern), ]
    AP <- data.frame[grep("A-P", data.frame$Pattern), ]
    write.table(data.frame, file=paste0(GeneName, "_AllPatternsGenes.csv"), quote=F, sep=",", row.names=F, col.names=T)
    write.table(Ectoderm, file=paste0(GeneName, "_EctodermGenes.csv"), quote=F, sep=",", row.names=F, col.names=T)
    write.table(Endoderm, file=paste0(GeneName, "_EndodermGenes.csv"), quote=F, sep=",", row.names=F, col.names=T)
    write.table(Mesoderm, file=paste0(GeneName, "_MesodermGenes.csv"), quote=F, sep=",", row.names=F, col.names=T)
    write.table(Ubiquitous, file=paste0(GeneName, "_UbiquitousGenes.csv"), quote=F, sep=",", row.names=F, col.names=T)

}

