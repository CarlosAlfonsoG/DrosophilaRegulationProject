###FindOverlapsinGenomicRangesKeepingmetadata##
## Cite: answerd by Valerie Obenchain Bioconductor: https://support.bioconductor.org/p/54470/
GrangesOverlapsWithMetadata <- function(Gr1, Gr2, GeneName){
          require(GenomicRanges); require(gsubfn); require(stringi)
          ranges <-subsetByOverlaps(Gr1, Gr2)
          hits <- findOverlaps(Gr1, Gr2)
          Pattern <- CharacterList(split(Gr2$Pattern[subjectHits(hits)], queryHits(hits)))
          ReporterIntensity <- CharacterList(split(Gr2$LevelofExp[subjectHits(hits)], queryHits(hits)))
          VTID <- CharacterList(split(Gr2$VTID[subjectHits(hits)], queryHits(hits)))
          mcols(ranges) <- DataFrame(mcols(ranges), Pattern, ReporterIntensity,VTID)
          AnnotatedGrangesWithPattern <<- ranges
          AnnotatedGrangesWithPattern
          #CleanData and generate Duplicate Rows for pattern
          AnnotatedGrangesWithPattern.df <<- as.data.frame(AnnotatedGrangesWithPattern)
          AnnotatedGrangesWithPattern.df$Pattern <<- as.character(AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <<- gsub("\\(", "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <<- gsub("\\)", "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <<- gsub('"', "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$Pattern <<- gsub("^c", "", AnnotatedGrangesWithPattern.df$Pattern)
          AnnotatedGrangesWithPattern.df$ID <<- seq.int(nrow(AnnotatedGrangesWithPattern.df))
          df <- AnnotatedGrangesWithPattern.df
          s <- strsplit(df$Pattern, split = ",")
          e <- data.frame(ID = rep(df$ID, sapply(s, length)), Pattern = unlist(s))
          colnames(e)[1] <- "ID"
          f <- merge(df, e, by="ID", all.y = T)
          f$ReporterIntensity <- as.character(f$ReporterIntensity)
          f$ReporterIntensity <<- stri_extract_first_regex(f$ReporterIntensity, "[0-9]")
          
          ##MakeGrangesToBed by Devon Ryan Biostars: https://www.biostars.org/p/89341/
          write.table(AnnotatedGrangesWithPattern, file=paste0(GeneName, "AnnotatedGrangesWithPattern.bed"), quote=F, sep="\t", row.names=F, col.names=F)
}

