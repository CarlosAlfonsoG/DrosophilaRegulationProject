library(GenomicRanges)
# generate 10 random segments
regions.gr <- GRanges(seqnames=Rle(rep(1,10)), ranges=IRanges(start=sample(1:1000,10), width=100), strand=Rle(rep("*",10)))

# generate 5 random genes
fiveGeneNames <- paste(sample(letters,5), sample(1:100,5), sep="")
genes.gr <- GRanges(seqnames=Rle(rep(1,5)), ranges=IRanges(start=sample(1:1000,5), width=100, names=fiveGeneNames), strand=Rle(sample(c("+","-"),5, replace=T)))

# what regions overlap what genes?
overlapGenes <- findOverlaps(regions.gr, genes.gr)

# Return any genes with an overlap.
# Convert the resulting "Hits" object to a data frame
# and use it as an index
overlapGenes.df <- as.data.frame(overlapGenes)
names(genes.gr[overlapGenes.df$subjectHits])

# extract the regions that hit genes
regionsWithHits <- as.data.frame(regions.gr[overlapGenes.df$queryHits])
# add the names of the genes
regionsWithHits$genes <- names(genes.gr)[overlapGenes.df$subjectHits]
#write output
write.table(regionsWithHitsKrup, file="Krupp_StarkMatch.bed", quote=F, sep="\t", row.names=F, col.names=F)