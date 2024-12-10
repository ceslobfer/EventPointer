sampleWhichCount <- function(sample_info, list_features) {
  list_which <- unlist(lapply(list_features, function(chrFeature) {
    if (length(chrFeature) != 0) {
      which <- range(chrFeature)
      seq <- as.character(unique(seqnames(which)))
      newWhich <- GRanges(
        seqnames = Rle(c(seq), c(1)),
        ranges = IRanges(min(start(which)), end = max(end(which))),
        strand = Rle(strand(c("-")), c(1))
      )
      newWhich
    }
  }))
  
  res <- c()
  for (sample in c(1:length(sample_info$file_bam))) {
    res <- c(
      res,
      lapply(
        X = list_which,
        FUN = listCHRFileCount,
        valSample = sample_info[sample, ],
        list_features = list_features
      )
    )
    
  }
  return(res)
}

listCHRFileCount <- function(which, valSample, list_features) {
  seq <- as.character(seqnames(which))
  features <- list_features[[seq]]
  file_bam <- valSample$file_bam
  paired_end <- valSample$paired_end
  sample_name <- valSample$sample_name
  addWhich <-
    list(file_bam, paired_end, sample_name, which, features)
  names(addWhich) <-
    c("file_bam", "paired_end", "sample_name", "which", "features")
  return(addWhich)
}


getSGFeatureCounts <- function(sample_info,
                               features,
                               min_anchor = 6,
                               counts_only = FALSE,
                               retain_coverage = FALSE,
                               verbose = FALSE,
                               cores = 1)
{
  
  if (!is(features, "SGFeatures")) {
    stop("features must be an SGFeatures object")
    
  }
  #Get which from features
  
  mergeFeatures <- features
  list_range <- range(mergeFeatures)
  hits <- findOverlaps(mergeFeatures, list_range)
  list_index <- split(queryHits(hits), subjectHits(hits))
  mergeFeatures <- mergeFeatures[queryHits(hits)]
  list_features <- split(mergeFeatures, seqnames(mergeFeatures))
  
  list_element_count <- sampleWhichCount(sample_info, list_features)
  clC <- makePSOCKcluster(cores, outfile = "reportCount.txt")
  list_counts <- clusterApplyLB(
    cl = clC,
    x = list_element_count,
    fun = getSGFeatureCountsTotal,
    min_anchor = min_anchor,
    retain_coverage = retain_coverage,
    verbose = verbose
  )
  list_counts <- unlist(list_counts, recursive = FALSE)
  
  on.exit(stopCluster(cl = clC))
  CountPerSample <- list()
  for (sample in sample_info$sample_name) {
    listCountSample <- list()
    for (chr in c(1:length(list_counts))) {
      nameCHR <- names(list_counts[chr])
      chr <- list_counts[[chr]]
      listStrand <- list()
      for (chrStrand in chr) {
        if (sample == chrStrand$sample_name) {
          listStrand <- c(listStrand, list(chrStrand$counts))
          
        }
      }
      if (length(listStrand) != 0) {
        listCountSample[[nameCHR]] <- listStrand
      }
      
    }
    
    listCountSample <-
      listCountSample[unique(as.character(seqnames(mergeFeatures)))]
    counts <- c()
    if (retain_coverage) {
      counts <- do.call(rbind, listCountSample)
      counts <- counts[order(unlist(list_index)),]
      
    } else {
      counts <- unlist(listCountSample, use.names = FALSE)
      counts <- counts[order(unlist(list_index))]
      
    }
    if (length(CountPerSample) == 0) {
      CountPerSample <- list(counts)
    } else{
      CountPerSample <- c(CountPerSample, list(counts))
    }
    
    
  }
  
  counts <- do.call(cbind, CountPerSample)
  
  if (counts_only)
    return(counts)
  
  sgfc <- makeSGFeatureCounts(
    rowRanges = features,
    colData = sample_info,
    counts = counts,
    min_anchor = min_anchor
  )
  
  return(sgfc)
  
}

getSGFeatureCountsTotal <- function(mainCount,min_anchor, retain_coverage, verbose)
{
  features <- mainCount$features
  sample_name <- mainCount$sample_name
  file_bam <- mainCount$file_bam
  paired_end <- mainCount$paired_end
  sample_name <- mainCount$sample_name
  which <- mainCount$which
  whichMod<- which
  chrName <- as.character(which@seqnames)
  seqlevel <- as.character(seqnames(which))
  strand <- unique(as.character(strand(features)))
  print(paste("READ CHR",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
  pairGap <- readGapPair(file_bam, paired_end, which, sample_name, verbose)
  listCount <- list()
  if (length(strand)==2) {
    featuresPlus <- features[strand(features) == "+"]
    featuresMinus <- features[strand(features) == "-"]
    print(paste("Start + strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    listCount[[chrName]] <- list("+" = processCounts(pairGap$gapPlus, featuresPlus, strand = "+",
                                                     sample_name, min_anchor, retain_coverage, verbose))
    print(paste("End + strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    
    gc()
    print(paste("Start - strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    listCount[[chrName]] <- c(listCount[[chrName]],list("-" = processCounts(pairGap$gapMinus, featuresMinus, strand = "-",
                                                                            sample_name, min_anchor, retain_coverage, verbose)))
    print(paste("End - strand",as.character(seqnames(which)), sample_name, Sys.time(), sep = " "))
    
    gc()
  }else{
    if (strand == "+"){
      featuresPlus <- features[strand(features) == "+"]
      listCount[[chrName]] <- list("+" = processCounts(pairGap$gapPlus, featuresPlus, strand = "+",
                                                       sample_name, min_anchor, retain_coverage, verbose))
      gc()
    }else{
      featuresMinus <- features[strand(features) == "-"]
      listCount[[chrName]] <- list("-" = processCounts(pairGap$gapMinus, featuresMinus, strand = "-",
                                                       sample_name, min_anchor, retain_coverage, verbose))
      gc()
    }
  }
  
  return(listCount)
  
}

processCounts <- function(gap, features, strand,
                          sample_name, min_anchor, retain_coverage, verbose){
  
  
  which <- range(features)
  chrName <- as.character(which@seqnames)
  if (isEmpty(gap$frag_exonic)) {
    res_count <- list(counts = as.integer(rep.int(0,length(features))), sample_name=sample_name)
    rm(gap)
    gc()
  }else{
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron
    
    rm(gap)
    gc()
    
    rangeExon <- unlist(range(frag_exonic))
    
    selectionRangeExon <- range(which(end(rangeExon) >= min(start(features)) & start(rangeExon) <= max(end(features))))
    print(c(selectionRangeExon, chrName, strand))
    if (is.infinite(selectionRangeExon[1]) | is.infinite(selectionRangeExon[2]) ) {
      rm(frag_exonic)
      rm(frag_intron)
      gc()
      res_count <- list(counts = as.integer(rep.int(0,length(features))), sample_name=sample_name)
    }else{
      frag_exonic <- frag_exonic[c(selectionRangeExon[1]:selectionRangeExon[2])]
      frag_intron <- frag_intron[c(selectionRangeExon[1]:selectionRangeExon[2])]
      
      
      
      ir <- SGSeq:::extractRangesFromFeatures(features)
      
      ## extract feature type and spliced boundaries
      type <- mcols(ir)$type
      spliceL <- mcols(ir)$spliceL
      spliceR <- mcols(ir)$spliceR
      
      i_J <- which(type == "J")
      i_E <- which(type == "E")
      i_S <- which(type %in% c("spliceL", "spliceR"))
      
      N <- rep(NA_integer_, length(ir))
      if (length(i_J) > 0) {
        
        N[i_J] <- junctionCompatible(ir[i_J], frag_exonic, frag_intron,
                                     min_anchor)
        
      }
      
      if (length(i_E) > 0) {
        
        E_index <- exonCompatible(ir[i_E], spliceL[i_E], spliceR[i_E],
                                  frag_exonic, frag_intron, FALSE)
        
        N[i_E] <- elementNROWS(E_index)
        
      }
      if (length(i_S) > 0) {
        
        N[i_S] <- splicesiteOverlap(ir[i_S],
                                    sub("splice", "", type[i_S], fixed = TRUE),
                                    frag_exonic, frag_intron, min_anchor, "unspliced")
        
      }
      if (retain_coverage) {
        
        counts <- DataFrame(N = N)
        counts$N_splicesite <- IntegerList(vector("list", nrow(counts)))
        counts$coverage <- RleList(IntegerList(vector("list", nrow(counts))))
        if (length(i_J) > 0) {
          
          counts$N_splicesite[i_J] <- splicesiteCounts(ir[i_J],
                                                       frag_exonic, frag_intron, min_anchor, "junction", "all")
          
        }
        if (length(i_E) > 0) {
          
          counts$N_splicesite[i_E] <- splicesiteCounts(ir[i_E],
                                                       frag_exonic, frag_intron, min_anchor, "exon", "spliced")
          counts$coverage[i_E] <- exonCoverage(ir[i_E], E_index,
                                               frag_exonic)
          
        }
        if (strand == "-") {
          
          counts$N_splicesite <- endoapply(counts$N_splicesite, rev)
          counts$coverage <- endoapply(counts$coverage, rev)
        }
        
      } else {
        
        counts <- N
        
      }
      if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))
      rm(N)
      rm(E_index)
      
      res_count <- list(counts = counts, sample_name=sample_name)
    }
  }
  return(res_count)
  
}


# rlang::env_unlock(env = asNamespace('SGSeq'))
# rlang::env_binding_unlock(env = asNamespace('SGSeq'))
# assign('getSGFeatureCounts', getSGFeatureCounts, envir = asNamespace('SGSeq'))
# assign('sampleWhichCount', sampleWhichCount, envir = asNamespace('SGSeq'))
# assign('listCHRFileCount', listCHRFileCount, envir = asNamespace('SGSeq'))
# assign('getSGFeatureCountsTotal', getSGFeatureCountsTotal, envir = asNamespace('SGSeq'))
# assign('processCounts', processCounts, envir = asNamespace('SGSeq'))
# rlang::env_binding_lock(env = asNamespace('SGSeq'))
# rlang::env_lock(asNamespace('SGSeq'))

