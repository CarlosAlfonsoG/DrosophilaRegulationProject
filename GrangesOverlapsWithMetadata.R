###FindOverlapsinGenomicRangesKeepingmetadata##
## Cite: answerd by Valerie Obenchain Bioconductor: https://support.bioconductor.org/p/54470/
GrangesOverlapsWithMetadata <- function(Gr1, Gr2, GeneName){
          require(GenomicRanges)
          ranges <-subsetByOverlaps(Gr1, Gr2)
          hits <- findOverlaps(Gr1, Gr2)
          Pattern <- CharacterList(split(Gr2$Pattern[subjectHits(hits)], queryHits(hits)))
          ReporterIntensity <- CharacterList(split(Gr2$LevelofExp[subjectHits(hits)], queryHits(hits)))
          VTID <- CharacterList(split(Gr2$VTID[subjectHits(hits)], queryHits(hits)))
          mcols(ranges) <- DataFrame(mcols(ranges), Pattern, ReporterIntensity,VTID)
          AnnotatedGrangesWithPattern <<- ranges
          AnnotatedGrangesWithPattern
          ##MakeGrangesToBed by Devon Ryan Biostars: https://www.biostars.org/p/89341/
          write.table(AnnotatedGrangesWithPattern, file=paste0(GeneName, "AnnotatedGrangesWithPattern.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}
