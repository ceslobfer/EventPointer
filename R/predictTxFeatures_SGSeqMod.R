predictTxFeatures <- function(sample_info,
                              which = NULL,
                              alpha = 2,
                              psi = 0,
                              beta = 0.2,
                              gamma = 0.2,
                              include_counts = FALSE,
                              retain_coverage = FALSE,
                              junctions_only = FALSE,
                              min_junction_count = NULL,
                              min_anchor = 1,
                              max_complexity = 20,
                              min_n_sample = 1,
                              min_overhang = NA,
                              verbose = FALSE,
                              cores = 1)
{
  SGSeq:::checkSampleInfo(sample_info)
  cat("Amo a ve")
  list_whichSamples <- c()
  if (is.null(which)) {
    which <- globalWhich(sample_info)
    list_whichSamples <-
      sampleWhichPredict(sample_info, alpha, min_junction_count, which, TRUE)
  } else {
    which <- expandUnstrandedRanges(which)
    list_whichSamples <-
      sampleWhichPredict(sample_info, alpha, min_junction_count, which, FALSE)
  }
  clF <- makeCluster(cores, outfile = "reportPrediction.txt")
  # clF <- makePSOCKcluster(cores, outfile = "reportPrediction.txt")
  # list_features <- clusterApplyLB(
  #   cl = clF,
  #   x = list_whichSamples,
  #   fun = predictTxFeaturesTotal,
  #   psi = psi,
  #   beta = beta,
  #   gamma = gamma,
  #   min_anchor = min_anchor,
  #   include_counts = include_counts,
  #   retain_coverage = retain_coverage,
  #   junctions_only = junctions_only,
  #   max_complexity = max_complexity,
  #   verbose = verbose
  # )
  list_features <- parLapplyLB(
    cl = clF,
    X = list_whichSamples,
    fun = predictTxFeaturesTotal,
    psi = psi,
    beta = beta,
    gamma = gamma,
    min_anchor = min_anchor,
    include_counts = include_counts,
    retain_coverage = retain_coverage,
    junctions_only = junctions_only,
    max_complexity = max_complexity,
    verbose = verbose
  )
  on.exit(stopCluster(cl = clF))
  featuresPerSample <- c()
  for (sample in sample_info$sample_name) {
    lFeaturesSample <- lapply(list_features, function(x) {
      if (sample == x$sample_name) {
        x$gr
      }
    })
    
    lFeaturesSample <- unlist(lFeaturesSample)
    lFeaturesSample <-
      lFeaturesSample[!vapply(lFeaturesSample, is.null, logical(1))]
    
    if (length(lFeaturesSample) == 0) {
      features <- TxFeatures()
      
    } else {
      features <- do.call(c, setNames(lFeaturesSample, NULL))
      features <- sort(features)
      features <- TxFeatures(features)
    }
    featuresPerSample <- c(featuresPerSample, features)
    
  }
  
  features <-
    mergeTxFeatures(featuresPerSample, min_n_sample = min_n_sample)
  
  if (!is.null(min_overhang)) {
    features <- processTerminalExons(features, min_overhang)
    
  }
  
  return(features)
  
}

expandUnstrandedRanges <- function(x)
{
  
  i <- which(strand(x) == "*")
  
  if (length(i) > 0) {
    
    additional <- x[i]
    strand(additional) <- "-"
    strand(x)[i] <- "+"
    x <- c(x, additional)
    
  }
  
  return(x)
  
}

validationParameters <-
  function(sample_info,
           alpha,
           min_junction_count,
           which) {
    file_bam = sample_info$file_bam
    paired_end = sample_info$paired_end
    read_length = sample_info$read_length
    frag_length = sample_info$frag_length
    lib_size = sample_info$lib_size
    sample_name = sample_info$sample_name
    
    if (is.null(min_junction_count) && is.null(alpha)) {
      stop("Need to provide min_junction_count or alpha")
      
    }
    
    if (is.null(min_junction_count) &&
        (is.null(read_length) ||
         is.null(frag_length) || is.null(lib_size))) {
      stop("For use with alpha, need to provide read_length,
            frag_length and lib_size")
      
    }
    
    if (!is.null(which) && !is(which, "GRanges")) {
      stop("argument which must be NULL or a GRanges object")
      
    }
    
    if (is.null(min_junction_count)) {
      min_junction_count <- SGSeq:::convertFpkmToCount(alpha, paired_end,
                                                       read_length, frag_length, lib_size)
      
    }
    
    return(
      list(
        min_junction_count = min_junction_count,
        file_bam = file_bam,
        paired_end = paired_end,
        sample_name = sample_name
      )
    )
  }


globalWhich <- function(sample_info) {
  res <- c()
  for (i in c(1:length(sample_info$file_bam))) {
    file_bam <- seqinfo(BamFile(sample_info$file_bam[i]))
    sl <- rep(seqlevels(file_bam), 2)
    st <- rep(c("+", "-"), rep(length(file_bam), 2))
    which <- GRanges(sl, IRanges(1, seqlengths(file_bam)[sl]), st)
    if (length(res) == 0) {
      res <- which
    } else{
      res <- c(res, which)
    }
  }
  res <- reduce(res)
  return(res)
  
}

sampleWhichPredict <-
  function(sample_info,
           alpha,
           min_junction_count,
           which,
           novo) {
    res <- c()
    for (sample in c(1:length(sample_info$file_bam))) {
      valSample <- validationParameters(sample_info[sample, ],
                                        alpha, min_junction_count, which)
      paired_end = valSample$paired_end
      searchWhich <- which
      if (paired_end & novo) {
        searchWhich <- searchWhich[searchWhich@strand == "+"]
      }
      bam_index <- idxstatsBam(file = sample_info[sample, ]$file_bam)
      res <-
        c(
          res,
          lapply(
            X = split(searchWhich, seq_along(searchWhich)),
            FUN = listCHRFilePredict,
            valSample = valSample,
            bam_index = bam_index
          )
        )
      
    }
    return(res)
  }

listCHRFilePredict <- function(range, valSample, bam_index) {
  file_bam = valSample$file_bam
  paired_end = valSample$paired_end
  sample_name = valSample$sample_name
  as.character(range@seqnames)
  index_chr <-
    bam_index[which(bam_index$seqnames == as.character(range@seqnames)), ]
  coefComplex <- index_chr$mapped / index_chr$seqlength
  min_junction_count = valSample$min_junction_count
  addWhich <-
    list(file_bam,
         paired_end,
         sample_name,
         min_junction_count,
         range,
         coefComplex)
  names(addWhich) <-
    c(
      "file_bam",
      "paired_end",
      "sample_name",
      "min_junction_count",
      "which",
      "coefComplex"
    )
  return(addWhich)
}

predictTxFeaturesTotal <- function(sWhich,file_bam, paired_end, which,
                                   min_junction_count, psi, beta, gamma, min_anchor, include_counts,
                                   retain_coverage, junctions_only, max_complexity, sample_name, verbose){
  which <- sWhich$which
  sample_name <- sWhich$sample_name
  seqlevel <- as.character(seqnames(which))
  
  file_bam <- sWhich$file_bam
  paired_end <- sWhich$paired_end
  which <- sWhich$which
  sample_name <- sWhich$sample_name
  min_junction_count <- sWhich$min_junction_count
  
  if (is(file_bam, "BamFile")) {
    
    si <- seqinfo(file_bam)
    
  } else {
    
    si <- seqinfo(BamFile(file_bam))
    
  }
  
  seqlevel <- as.character(seqnames(which))
  strand <- as.character(strand(which))
  #print(paste(c(seqlevel, strand), sep = " "))
  
  if (paired_end) {
    pairGap <- readGapPair(file_bam, paired_end, which, sample_name, verbose)
    
    print(paste(sample_name,"read CHR ", seqlevel, Sys.time(), sep = " "))
    
    gc(full = TRUE)
    
    frag_exonic <- pairGap$gapPlus$frag_exonic
    frag_intron <- pairGap$gapPlus$frag_intron
    grPlus <- constructPaired(frag_exonic, frag_intron, min_junction_count,
                              psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                              junctions_only, max_complexity, sample_name, seqlevel, "+", si)
    print(paste("    ",sample_name,seqlevel, "+",Sys.time(), sep = " "))
    
    frag_exonic <- pairGap$gapMinus$frag_exonic
    frag_intron <- pairGap$gapMinus$frag_intron
    grMinus <- constructPaired(frag_exonic, frag_intron, min_junction_count,
                               psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                               junctions_only, max_complexity, sample_name, seqlevel, "-", si)
    
    print(paste("    ",sample_name,seqlevel, "-",Sys.time(), sep = " "))
    rm(frag_exonic)
    rm(frag_intron)
    rm(pairGap)
    gc()
    gr <- c(grPlus,grMinus)
  }else{
    gap <- readGAlignments(file_bam, paired_end, which, sample_name, verbose)
    frag_exonic <- gap$frag_exonic
    frag_intron <- gap$frag_intron
    gr <- constructPaired(frag_exonic, frag_intron, min_junction_count,
                          psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                          junctions_only, max_complexity, sample_name, seqlevel, strand, si)
    rm(frag_exonic)
    rm(frag_intron)
    rm(gap)
  }
  
  
  gc(full=TRUE)
  if (verbose) generateCompleteMessage(paste(sample_name, gr2co(which)))
  res <- list(sample_name=sample_name,gr=gr)
  
  return(res)
  
}
predictJunctions <- function(frag_exonic, frag_intron,
                             min_junction_count, psi, min_anchor, retain_coverage)
{
  
  ## extract all splice junctions
  junctions <- unique(unlist(frag_intron)) + 1
  if (length(junctions) == 0) { return() }
  mcols(junctions) <- DataFrame("type" = rep("J", length(junctions)))
  ## consider splice junctions with counts at least min_junction_count
  #print("Start junctinoCompatible")
  mcols(junctions)$N <- junctionCompatible(junctions, frag_exonic,
                                           frag_intron, min_anchor)
  
  junctions <- junctions[which(mcols(junctions)$N >= min_junction_count)]
  
  if (length(junctions) == 0) { return() }
  
  ## consider splice junctions and splicesites with
  ## counts >= psi * max(splicesite counts)
  
  ## Note left/right (L/R) nomenclature: for the LHS splice site,
  ## the spliced boundary is situated on the right, for the RHS
  ## splice site, the spliced boundary is situated on the left
  
  if (psi > 0 || retain_coverage) {
    mcols(junctions)$N_splicesite <- splicesiteCounts(junctions,
                                                      frag_exonic, frag_intron, min_anchor, "junction", "all")
    junctions <- junctions[which(mcols(junctions)$N >=
                                   psi * max(mcols(junctions)$N_splicesite))]
    if (length(junctions) == 0) { return() }
    
  }
  
  if (!retain_coverage) {
    
    mcols(junctions)$N_splicesite <- NULL
    
  }
  junctions <- SGSeq:::completeMcols(junctions, retain_coverage)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(junctions)
  
}

readGapPair <- function(file, paired_end, which, sample_name, verbose)
{
  
  if (length(which) != 1) {
    
    stop("which argument must have length 1")
    
  }
  
  flag <- scanBamFlag(isSecondaryAlignment = FALSE)
  param <- ScanBamParam(flag = flag, tag = "XS", which = which)
  if (paired_end) {
    gap <- suppressWarnings(readGAlignmentPairs(file = file,
                                                param = param)) 
    gap <- SGSeq:::propagateXS(gap)
  }else{
    gap <- suppressWarnings(readGAlignments(file = file, param = param))
  }
  
  
  
  gap <- SGSeq:::filterGap(gap)
  
  mcols(gap)$strand <- SGSeq:::XS2strand(mcols(gap)$XS)
  
  gapPlus <- fragExonIntron(gap = gap, strand="+",verbose)
  gapMinus <- fragExonIntron(gap = gap, strand="-",verbose)
  rm(gap)
  gc()
  return(list(gapPlus=gapPlus,gapMinus=gapMinus))
  
}
generateWarningMessage <- function (fun_name, item, msg) 
{
  message(makeWarningMessage(fun_name, item, msg))
}
fragExonIntron <- function(gap,strand,verbose){
  gap <- gap[mcols(gap)$strand %in% c(strand, "*")]
  frag_exonic <- reduce(ranges(grglist(gap, drop.D.ranges = TRUE)))
  frag_intron <- ranges(junctions(gap))
  diff <- setdiff(frag_exonic, frag_intron)
  excl <- which(sum(width(frag_exonic)) > sum(width(diff)))
  
  if (length(excl) > 0) {
    
    if (verbose) {
      
      msg <- paste(
        "filtered",
        length(excl),
        "inconsistent paired alignments in",
        gr2co(which))
      
      SGSeq:::generateWarningMessage(
        "readGap",
        sample_name,
        msg)
      
    }
    
    frag_exonic <- frag_exonic[-excl]
    frag_intron <- frag_intron[-excl]
    
  }
  
  gap <- list(frag_exonic = frag_exonic, frag_intron = frag_intron)
  rm(frag_exonic)
  rm(frag_intron)
  return(gap)
}

togroup0 <- S4Vectors:::quick_togroup

constructPaired <- function(frag_exonic, frag_intron, min_junction_count,
                            psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                            junctions_only, max_complexity, sample_name, seqlevel, strand, si){
  
  if (length(frag_exonic) == 0) {
    
    gr <- NULL
    
  } else {
    library(profvis)
    ir <- predictSpliced(frag_exonic, frag_intron, min_junction_count,
                         psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                         junctions_only, max_complexity, sample_name, seqlevel, strand)
    
    if (is.null(ir)) {
      gr <- NULL
      
    } else {
      gr <- constructGRangesFromRanges(ir, seqlevel, strand, si)
    }
    
  }
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(gr)
}

##' Ranges-based identification of splice junctions and exons.
##'
##' @title Ranges-based identification of splice junctions and exons
##' @inheritParams predictTxFeaturesPerSample
##' @param frag_exonic \code{IRangesList} with exonic regions from alignments
##' @param frag_intron \code{IRangesList} with introns implied by spliced
##'   alignments
##' @param min_anchor Integer specifiying minimum anchor length
##' @param seqlevel \code{seqlevel} to be processed
##' @param strand \code{strand} to be processed
##' @return \code{IRanges} with predicted features
##' @keywords internal
##' @author Leonard Goldstein

co2str <- function (seqlevel, start, end, strand) 
{
  paste0(seqlevel, ":", start, "-", end, ":", strand)
}

predictSpliced <- function(frag_exonic, frag_intron, min_junction_count,
                           psi, beta, gamma, min_anchor, include_counts, retain_coverage,
                           junctions_only, max_complexity, sample_name, seqlevel, strand)
{
  junctions <- predictJunctions(frag_exonic, frag_intron,
                                min_junction_count, psi, min_anchor, retain_coverage)
  if (is.null(junctions)) { return() }
  
  if (!junctions_only) {
    lower <- max(min_junction_count * beta, 1)
    frag_coverage <- coverage(unlist(frag_exonic))
    islands <- slice(frag_coverage, lower, rangesOnly = TRUE)
    ## skip problematic regions
    
    if (!is.na(max_complexity)) {
      
      ir <- as(slice(coverage(junctions), max_complexity), "IRanges")
      
      if (length(ir) > 0) {
        
        junctions_stripped <- junctions
        mcols(junctions_stripped) <- NULL
        loci <- reduce(c(junctions_stripped, islands))
        excl <- loci[loci %over% ir]
        
        junctions <- junctions[!junctions %over% excl]
        islands <- islands[!islands %over% excl]
        
        excl_str <- co2str(seqlevel, start(excl), end(excl), strand)
        
        SGSeq:::generateWarningMessage(
          "predictSpliced",
          sample_name,
          paste("skipping", excl_str))
        
        if (length(junctions) == 0) { return() }
        
      }
      
    }
    
    features <- junctions
    # print("Start extractSplicesitesFromJunctions")
    splicesites_L <- extractSplicesitesFromJunctions(junctions, "L")
    splicesites_R <- extractSplicesitesFromJunctions(junctions, "R")
    splicesites <- c(splicesites_L, splicesites_R)
    # print("Start predictCandidatesInternal")
    candidates <- predictCandidatesInternal(islands, splicesites,
                                            frag_coverage, beta)
    # print("Start predictExonsInternal")
    exons_I <- predictExonsInternal(candidates, frag_exonic,
                                    frag_intron, beta, min_anchor, include_counts, retain_coverage)
    if (!is.null(exons_I)) { features <- c(features, exons_I) }
    # print("Start predictCandidatesTerminal L")
    candidates <- predictCandidatesTerminal(islands, splicesites, "exon_L")
    exons_L <- predictExonsTerminal(candidates, frag_exonic, frag_intron,
                                    gamma, min_anchor, "exon_L", include_counts, retain_coverage)
    if (!is.null(exons_L)) { features <- c(features, exons_L) }
    # print("Start predictCandidatesTerminal R")
    candidates <- predictCandidatesTerminal(islands, splicesites, "exon_R")
    exons_R <- predictExonsTerminal(candidates, frag_exonic, frag_intron,
                                    gamma, min_anchor, "exon_R", include_counts, retain_coverage)
    if (!is.null(exons_R)) { features <- c(features, exons_R) }
    
  } else {
    
    features <- junctions
    
  }
  
  if (!include_counts) {
    
    mcols(features)$N <- NULL
    
  }
  rm(frag_exonic)
  rm(frag_intron)
  rm(junctions)
  gc(full=TRUE)
  features <- sort(features)
  
  return(features)
  
}

extractSplicesitesFromJunctions <- function(junctions, type = c("L", "R"))
{
  
  type <- match.arg(type)
  S <- flank(junctions, -1, switch(type, "L" = TRUE, "R" = FALSE))
  S_pos <- as.character(start(S))
  pos_N <- tapply(mcols(junctions)$N, S_pos, sum)
  pos_N <- setNames(as.integer(pos_N), names(pos_N))
  i <- which(!duplicated(S_pos))
  S <- S[i]
  S_pos <- S_pos[i]
  mcols(S) <- DataFrame(type = rep(type, length(S)),
                        N = pos_N[match(S_pos, names(pos_N))])
  return(S)
  
}

predictExonsTerminal <- function(candidates, frag_exonic, frag_intron, relCov,
                                 min_anchor, type = c("exon_L", "exon_R"), include_counts, retain_coverage)
{
  
  type <- match.arg(type)
  
  if (length(candidates) == 0) { return() }
  
  spliceL <- switch(type, "exon_L" = FALSE, "exon_R" = TRUE)
  spliceR <- switch(type, "exon_L" = TRUE, "exon_R" = FALSE)
  
  index <- exonCompatible(candidates, spliceL, spliceR,
                          frag_exonic, frag_intron, FALSE)
  coverage <- exonCoverage(candidates, index, frag_exonic)
  
  splicesite <- flank(candidates, -1,
                      start = switch(type, "exon_L" = FALSE, "exon_R" = TRUE))
  N_splicesite <- splicesiteOverlap(splicesite,
                                    switch(type, "exon_L" = "R", "exon_R" = "L"),
                                    frag_exonic, frag_intron, min_anchor, "spliced")
  
  ranges <- (coverage >= relCov * N_splicesite)
  el <- elementNROWS(ranges)
  rl <- runLength(ranges)
  
  if (type == "exon_L") {
    
    ir <- IRanges(end = el, width = SGSeq:::plast(rl))
    
  }
  if (type == "exon_R") {
    
    ir <- IRanges(start = 1, width = SGSeq:::pfirst(rl))
    
  }
  
  if (length(ir) ==    0) { return() }
  
  exons <- SummarizedExperiment::shift(ir, start(candidates) - 1)
  mcols(exons) <- DataFrame(type = rep(type, length(exons)))
  
  if (include_counts) {
    
    mcols(exons)$N <- exonCompatible(exons, spliceL, spliceR,
                                     frag_exonic, frag_intron)
    
  }
  
  if (retain_coverage) {
    
    mcols(exons)$N_splicesite <- as(N_splicesite, "IntegerList")
    mcols(exons)$coverage <- coverage[
      setNames(split(ir, seq_along(ir)), NULL)]
    
  }
  
  exons <- SGSeq:::completeMcols(exons, retain_coverage)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(exons)
  
}
predictCandidatesTerminal <- function(islands, splicesites,
                                      type = c("exon_L", "exon_R"))
{
  
  type <- match.arg(type)
  
  splicesites <- splicesites[mcols(splicesites)$type ==
                               switch(type, "exon_L" = "L", "exon_R" = "R")]
  hits <- findOverlaps(splicesites, islands)
  spliced_boundary <- splicesites[queryHits(hits)]
  island <- islands[subjectHits(hits)]
  
  if (type == "exon_L") {
    
    candidates <- IRanges(start(island), end(spliced_boundary))
    
  }
  if (type == "exon_R") {
    
    candidates <- IRanges(start(spliced_boundary), end(island))
    
  }
  
  mcols(candidates) <- DataFrame(N = mcols(splicesites)$N[queryHits(hits)])
  
  return(candidates)
  
}
predictExonsInternal <- function(candidates, frag_exonic, frag_intron, relCov,
                                 min_anchor, include_counts, retain_coverage)
{
  
  if (length(candidates) == 0) { return() }
  # print("predictExonsInternal 1")
  candidate_index <- exonCompatible(candidates, TRUE, TRUE,
                                    frag_exonic, frag_intron, FALSE)
  # print("predictExonsInternal 2")
  candidate_coverage <- exonCoverage(candidates, candidate_index,
                                     frag_exonic)
  # print("predictExonsInternal 3")
  candidate_N_splicesite <- splicesiteCounts(candidates, frag_exonic,
                                             frag_intron, min_anchor, "exon", "spliced")
  # print("predictExonsInternal 4")
  index <- which(min(candidate_coverage) >=
                   relCov * min(candidate_N_splicesite))
  if (length(index) == 0) { return() }
  # print("predictExonsInternal 5")
  exons <- candidates[index]
  mcols(exons) <- DataFrame("type" = rep("I", length(exons)))
  
  if (include_counts) {
    
    mcols(exons)$N <- exonCompatible(exons, TRUE, TRUE, frag_exonic,
                                     frag_intron)
    
  }
  # print("predictExonsInternal 6")
  if (retain_coverage) {
    
    mcols(exons)$N_splicesite <- candidate_N_splicesite[index]
    mcols(exons)$coverage <- candidate_coverage[index]
    
  }
  
  exons <- SGSeq:::completeMcols(exons, retain_coverage)
  # print("predictExonsInternal 7")
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(exons)
  
}
predictCandidatesInternal <- function(islands, splicesites, frag_coverage,
                                      relCov)
{
  
  ## for each island, identify overlapping splice sites
  
  island_splicesite <- as.list(findOverlaps(islands, splicesites))
  
  ## for each island, obtain all pairs of overlapping splice sites
  
  island_splicesite_pairs <- mapply(expand.grid,
                                    island_splicesite, island_splicesite, SIMPLIFY = FALSE)
  # print("Aski 1")
  splicesite_pairs <- unique(do.call(SGSeq:::rbindDfsWithoutRowNames,
                                     island_splicesite_pairs))
  
  ## retain pairs of splice sites that are consistent with
  ## flanking an internal exon
  
  N_1 <- mcols(splicesites)$N[splicesite_pairs[, 1]]
  N_2 <- mcols(splicesites)$N[splicesite_pairs[, 2]]
  type_1 <- mcols(splicesites)$type[splicesite_pairs[, 1]]
  type_2 <- mcols(splicesites)$type[splicesite_pairs[, 2]]
  pos_1 <- start(splicesites)[splicesite_pairs[, 1]]
  pos_2 <- start(splicesites)[splicesite_pairs[, 2]]
  i <- which(type_1 == "R" & type_2 == "L" & pos_1 <= pos_2)
  candidates <- IRanges(pos_1[i], pos_2[i])
  mcols(candidates) <- DataFrame(N = IntegerList(
    mapply(c, N_1[i], N_2[i], SIMPLIFY = FALSE)))
  
  if (length(candidates) > 0) {
    
    ## retain candidate internal exons with sufficient read coverage
    
    candidates_frag_coverage <- split(frag_coverage[candidates],
                                      togroup0(candidates))
    i <- which(min(candidates_frag_coverage) >=
                 relCov * min(mcols(candidates)$N))
    candidates <- candidates[i]
    
  }
  gc(full=TRUE)
  return(candidates)
  
}

constructGRangesFromRanges <- function(x, seqname, strand, seqinfo)
{
  
  x_mcols <- mcols(x)
  mcols(x) <- NULL
  
  if (strand == "+") {
    
    x_mcols_type <- as.character(x_mcols$type)
    x_mcols_type <- sub("exon_L", "F", x_mcols_type, fixed = TRUE)
    x_mcols_type <- sub("exon_R", "L", x_mcols_type, fixed = TRUE)
    x_mcols$type <- factor(x_mcols_type,
                           levels = c("J", "I", "F", "L", "U"))
    
  }
  if (strand == "-") {
    
    x_mcols_type <- as.character(x_mcols$type)
    x_mcols_type <- sub("exon_L", "L", x_mcols_type, fixed = TRUE)
    x_mcols_type <- sub("exon_R", "F", x_mcols_type, fixed = TRUE)
    x_mcols$type <- factor(x_mcols_type,
                           levels = c("J", "I", "F", "L", "U"))
    
    if ("N_splicesite" %in% names(x_mcols)) {
      
      x_mcols$N_splicesite <- endoapply(x_mcols$N_splicesite, rev)
      
    }
    if ("coverage" %in% names(x_mcols)) {
      
      x_mcols$coverage <- endoapply(x_mcols$coverage, rev)
      
    }
    
  }
  
  gr <- GRanges(rep(seqname, length(x)), x, rep(strand, length(x)),
                x_mcols, seqinfo = seqinfo)
  
  rm(x)
  gc(full=TRUE)
  
  return(gr)
  
}

# rlang::env_unlock(env = asNamespace('SGSeq'))
# rlang::env_binding_unlock(env = asNamespace('SGSeq'))
# assign('predictTxFeatures', predictTxFeatures, envir = asNamespace('SGSeq'))
# assign('validationParameters', validationParameters, envir = asNamespace('SGSeq'))
# assign('globalWhich', globalWhich, envir = asNamespace('SGSeq'))
# assign('sampleWhichPredict', sampleWhichPredict, envir = asNamespace('SGSeq'))
# assign('listCHRFilePredict', listCHRFilePredict, envir = asNamespace('SGSeq'))
# assign('constructPaired', constructPaired, envir = asNamespace('SGSeq'))
# assign('predictSpliced', predictSpliced, envir = asNamespace('SGSeq'))
# assign('constructGRangesFromRanges', constructGRangesFromRanges, envir = asNamespace('SGSeq'))
# assign('predictTxFeaturesTotal', predictTxFeaturesTotal, envir = asNamespace('SGSeq'))
# rlang::env_binding_lock(env = asNamespace('SGSeq'))
# rlang::env_lock(asNamespace('SGSeq'))
