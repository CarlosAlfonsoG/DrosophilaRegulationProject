###FindOverlapsinGenomicRangesKeepingmetadata##
## Cite: answerd by Valerie Obenchain Bioconductor: https://support.bioconductor.org/p/54470/
GrangesOverlapsWithMetadata <- function(Gr1, Gr2, GeneName){
          require(GenomicRanges); require(gsubfn)
          ranges <-subsetByOverlaps(Gr1, Gr2)
          hits <- findOverlaps(Gr1, Gr2)
          Pattern <- CharacterList(split(Gr2$Pattern[subjectHits(hits)], queryHits(hits)))
          ReporterIntensity <- CharacterList(split(Gr2$LevelofExp[subjectHits(hits)], queryHits(hits)))
          VTID <- CharacterList(split(Gr2$VTID[subjectHits(hits)], queryHits(hits)))
          mcols(ranges) <- DataFrame(mcols(ranges), Pattern, ReporterIntensity,VTID)
          AnnotatedGrangesWithPattern <<- ranges
          AnnotatedGrangesWithPattern
          #CleanData and generate Duplicate Rows for pattern
          AnnotatedGrangesWithPattern.df <- as.data.frame(AnnotatedGrangesWithPattern)
          AnnotatedGrangesWithPattern.df$Pattern <- as.character(AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <- gsub("\\(", "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <- gsub("\\)", "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <- gsub('"', "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <- gsub("^c", "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$ID <- seq.int(nrow(AnnotatedGrangesWithPattern.df))
          s <- strsplit(df$V2, split = ",")
          e <- data.frame(ID = rep(A1, sapply(s, length)), Pattern = unlist(s))
          f <- merge(A1, e[, c("ID", "Pattern")], by="ID", all.y = T)
          
          ##MakeGrangesToBed by Devon Ryan Biostars: https://www.biostars.org/p/89341/
          write.table(AnnotatedGrangesWithPattern, file=paste0(GeneName, "AnnotatedGrangesWithPattern.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}

