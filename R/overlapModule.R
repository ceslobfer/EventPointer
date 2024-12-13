##' Identify fragments compatible with splice junctions.
##'
##' @title Compatible fragment counts for splice junctions
##' @inheritParams exonCompatible
##' @param junctions \code{IRanges} of splice junctions
##' @param min_anchor Integer specifying minimum anchor length
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

junctionCompatible <- function(junctions, frag_exonic, frag_intron,
                               min_anchor, counts = TRUE)
{
  
  frag_intron <- filterIntrons(frag_intron, frag_exonic, min_anchor)
  
  introns <- junctions - 1
  
  hits <- findOverlapsRanges(introns, frag_intron, "equal")
  junction_index <- as.list(hits)
  rm(frag_exonic)
  rm(frag_intron)
  rm(junctions)
  gc(full=TRUE)
  if (counts) elementNROWS(junction_index)
  else junction_index
  
}

filterIntrons <- function(frag_intron, frag_exonic, min_anchor)
{
  unlisted_intron <- unlist(frag_intron)
  
  f_L <- as(flank(unlisted_intron, min_anchor, TRUE), "IRangesList")
  f_L <- intersect(f_L, frag_exonic[togroup0(frag_intron)])
  w_L <- sum(width(f_L))
  f_R <- as(flank(unlisted_intron, min_anchor, FALSE), "IRangesList")
  f_R <- intersect(f_R, frag_exonic[togroup0(frag_intron)])
  w_R <- sum(width(f_R))
  i <- which(w_L == min_anchor & w_R == min_anchor)
  
  filtered <- setNames(split(unlisted_intron[i],factor(togroup0(frag_intron)[i], seq_along(frag_intron))), NULL)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  return(filtered)
  
}

##' Identify fragments compatible with exons.
##'
##' @title Compatible fragment counts for exons
##' @param exons \code{IRanges} of exons
##' @param spliceL Logical vector indicating whether LHS boundary is spliced
##' @param spliceR Logical vector indicating whether RHS boundary is spliced
##' @param frag_exonic \code{IRangesList} of exonic regions, one entry
##'   per fragment
##' @param frag_intron \code{IRangesList} of introns, one entry per fragment
##' @param counts Logical indicating whether counts or indices of
##'   compatible fragments should be returned
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

dimMatrixFix <- function(matrix1, matrix2){
  if(dim(matrix1)[1] != dim(matrix2)[1]){
    if(nrow(matrix1) > nrow(matrix2)){
      rows_to_add <- nrow(matrix1) - nrow(matrix2)
      matrix2 <- rbind(matrix2,
                       sparseMatrix(i = NULL,j=NULL,dims = c(rows_to_add,ncol(matrix2)))
      )
    }
  }
  return(matrix2)
}

##' Identify fragments with alignments extending across exon/intron boundaries.
##'
##' @title Compatible fragment counts for splice sites
##' @inheritParams exonCompatible
##' @param splicesites \code{IRanges} of splice sites
##' @param side Character vector indicating whether the spliced boundary
##'   is to the left (\dQuote{L}) or right (\dQuote{R}) of the splice site
##' @param min_anchor Integer specifiying minimum anchor length
##' @param include Character string indicating whether considered fragments
##'   should be all that overlap the splice site (\dQuote{all}), those
##'   that are spliced at the site (\dQuote{spliced}) or those that are
##'   not spliced, i.e. extend into the adjacent intron (\dQuote{unspliced})
##' @return Counts or list of indices of compatible fragments
##' @keywords internal
##' @author Leonard Goldstein

splicesiteOverlap <- function(splicesites, side, frag_exonic, frag_intron,min_anchor, include = c("all", "spliced", "unspliced"), counts = TRUE){
  
  include <- match.arg(include)
  
  if (length(side) == 1) side <- rep(side, length(splicesites))
  
  start <- setNames(c(L = TRUE, R = FALSE)[side], NULL)
  flanking_intron <- flank(splicesites, 1, start = start)
  hits <- findOverlapsRanges(splicesites, frag_exonic, out="hits")
  
  if (include == "spliced" || include == "all") {
    
    frag_intron <- filterIntrons(frag_intron, frag_exonic, min_anchor)
    tmp <- findOverlapsRanges(flanking_intron, frag_intron, out="hits")
    hits_spliced <- intersect(hits, tmp)
    
  }
  
  if (include == "unspliced" || include == "all") {
    
    tmp <- findOverlapsRanges(flanking_intron, frag_exonic, out="hits")
    hits_unspliced <- intersect(hits, tmp)
    rm(tmp)
    qH <- queryHits(hits_unspliced)
    sH <- subjectHits(hits_unspliced)
    intron_anchor <- as(flank(splicesites, min_anchor, start = start),"IRangesList")
    exonic_anchor <- as(flank(splicesites, -min_anchor, start = start),"IRangesList")
    w_E <- sum(width(intersect(exonic_anchor[qH], frag_exonic[sH])))
    w_I <- sum(width(intersect(intron_anchor[qH], frag_exonic[sH])))
    hits_unspliced <- hits_unspliced[w_E == min_anchor & w_I == min_anchor]
    
  }
  
  if (include == "spliced") {
    hits <- hits_spliced
    rm(hits_spliced)
  } else if (include == "unspliced") {
    hits <- hits_unspliced
    rm(hits_unspliced)
    rm()
  } else if (include == "all") {
    hits <- union(hits_spliced, hits_unspliced)
    rm(hits_spliced)
    rm(hits_unspliced)
  }
  
  splicesite_index <- as.list(hits)
  rm(hits)
  rm(frag_exonic)
  rm(frag_intron)
  gc(full=TRUE)
  
  if (counts) elementNROWS(splicesite_index)
  else splicesite_index
}

##' Modified \code{findOverlaps} function for \code{IRanges},
##' \code{IRangesList} objects that behaves analogous to
##' \code{findOverlaps} for \code{GRanges}, \code{GRangesList} objects.
##'
##' @title Modified \code{findOverlaps} function for \code{IRanges},
##'   \code{IRangesList} objects
##' @param query \code{IRanges} or \code{IRangesList} object
##' @param subject \code{IRanges} or \code{IRangesList} object
##' @param type Passed to \code{findOverlaps}
##' @return \code{Hits} object
##' @keywords internal
##' @author Leonard Goldstein
##' 
##' 

exonCompatible<- function(exons, spliceL, spliceR, frag_exonic,frag_intron, counts = TRUE){
  
  if (length(spliceL) == 1) spliceL <- rep(spliceL, length(exons))
  if (length(spliceR) == 1) spliceR <- rep(spliceR, length(exons))
  res<- c()
  
  lenSubject <- length(frag_exonic)
  lenQuery <- length(exons)
  
  subject_unlisted <- unlist(frag_exonic)
  subject_togroup <- togroup0(frag_exonic)
  query_unlisted <- exons
  query_togroup <- seq_along(exons)
  
  df <- cbind(start(subject_unlisted), end(subject_unlisted))
  options(digits=12)
  uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2),nrow = 2))
 
  posiciones_repetidas <- split(seq_along(uniqueVector), uniqueVector)[as.character(unique(uniqueVector))]
  
  new_subject_unlist <- df[!duplicated(uniqueVector),]
  rm(df)
  rm(uniqueVector)
  if (length(nrow(new_subject_unlist)) == 0){
    if (length(new_subject_unlist) == 2) {
      new_subject_unlist <- IRanges(start = new_subject_unlist[1],end = new_subject_unlist[2])
    }else{
      new_subject_unlist <- IRanges(start = NULL,end = NULL)
    }
  }else{
    new_subject_unlist <- IRanges(start = new_subject_unlist[,1],end = new_subject_unlist[,2])
  }
  
  hits_exonic <- list()
  positionsVector <- split(c(1:lenQuery), ceiling(seq_along(query_unlisted)/3000))
  for(val in c(1:length(positionsVector))){
    pos <- positionsVector[val]
    hits_listed <- as.list(findOverlaps(query_unlisted[unlist(pos)], new_subject_unlist, type = "any"))
    if (length(hits_listed) == 0) {
      hits_exonic <- hits_listed
      rm(hits_listed)
    }else{
      hits_exonic <- c(hits_exonic,hits_listed)
      rm(hits_listed)
    }
    if(val%%10 == 0){
      gc()
    }
  }
  
  gc()
  
  hits_introns <- findOverlapsRanges(exons, frag_intron)
  
  excl1 <- findOverlapsRanges(flank(exons, 1, TRUE), frag_exonic)
  excl1[which(!spliceL)] <- 0L
  
  excl2 <- findOverlapsRanges(flank(exons, 1, FALSE), frag_exonic)
  excl2[which(!spliceR)] <- 0L
  lastErase <- 1
  if (counts){
    res <- lapply(seq_along(excl1),function(i){
      eraseGroup <- c(hits_introns[[i]],excl1[[i]],excl2[[i]])
      lgroup <- length(setdiff(subject_togroup[unlist(posiciones_repetidas[hits_exonic[[i]]])],eraseGroup))
      if(i%%10000 == 0){
        excl1[c(lastErase:i)] <<- 0L
        excl2[c(lastErase:i)] <<- 0L
        hits_introns[c(lastErase:i)] <<- 0L
        lastErase <- i
        gc()
      }
      
      return(lgroup)
    })
    
  }else{
    res <- lapply(seq_along(excl1),function(i){
      eraseGroup <- c(hits_introns[[i]],excl1[[i]],excl2[[i]])
      lgroup <- setdiff(subject_togroup[unlist(posiciones_repetidas[hits_exonic[[i]]])],eraseGroup)
      res <- c(res, list(lgroup))
      if(i != 1 & i%%10000 == 0){
        excl1[c(lastErase:i)] <<- 0L
        excl2[c(lastErase:i)] <<- 0L
        hits_introns[c(lastErase:i)] <<- 0L
        lastErase <<- i
        gc()
      }
      return(lgroup)
    })
    
  }
  
  rm(excl1)
  rm(excl2)
  rm(hits_introns)
  gc()
  
  rm(new_subject_unlist)
  rm(subject_togroup)
  rm(query_unlisted)
  rm(subject_unlisted)
  rm(posiciones_repetidas)
  
  gc()
  
  return(res)
}

findOverlapsRanges <- function(query, subject, type = "any", out = "list")
{
  
  subject_hits <- c()
  query_hits <- c()
  lenSubject <- length(subject)
  lenQuery <- length(query)
  
  if (is(query, "IRangesList")) {
    query_unlisted <- unlist(query)
    query_togroup <- togroup0(query)
  }else {
    query_unlisted <- query
    query_togroup <- seq_along(query)
  }
  if (is(subject, "IRangesList")) {
    subject_unlisted <- unlist(subject)
    subject_togroup <- togroup0(subject)
    rm(subject)
  }else {
    subject_unlisted <- subject
    subject_togroup <- seq_along(subject)
    rm(subject)
  }
  hits <- c()
  gc()
  if (type == "equal") {
    hits_unlisted <- findMatches(query_unlisted, subject_unlisted)
    subject_hits <- subjectHits(hits_unlisted)
    query_hits <- queryHits(hits_unlisted)
    qH <- query_togroup[query_hits]
    sH <- subject_togroup[subject_hits]
    df <- cbind(qH, sH)
    uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2),nrow = 2))
    df <- as.data.frame(df)
    hits <- df[!duplicated(uniqueVector),]
    
    rm(uniqueVector)
    rm(df)
    rm(sH)
    rm(qH)
    rm(query_hits)
    rm(subject_hits)
    
    rm(hits_unlisted)
    gc(full=TRUE)
    
    hits <- Hits(as.integer(hits[, 1]), as.integer(hits[, 2]),lenQuery, lenSubject, sort.by.query = TRUE)
    
  }else {
    if (!is(query, "IRangesList")) {
      hits <- modFindOverlap(query_unlisted,subject_unlisted,subject_togroup,lenQuery, lenSubject, out)
    }else{
      hits_unlisted <- findOverlaps(query_unlisted, subject_unlisted, type = "any")
      subject_hits <- subjectHits(hits_unlisted)
      query_hits <- queryHits(hits_unlisted)
      qH <- query_togroup[query_hits]
      sH <- subject_togroup[subject_hits]
      df <- cbind(qH, sH)
      uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2)))
      df <- as.data.frame(df)
      hits <- df[!duplicated(uniqueVector),]
      rm(uniqueVector)
      rm(df)
      rm(sH)
      rm(qH)
      rm(query_hits)
      rm(subject_hits)
      
      rm(hits_unlisted)
      gc(full=TRUE)
      
      hits <- Hits(as.integer(hits[, 1]), as.integer(hits[, 2]),lenQuery, lenSubject, sort.by.query = TRUE)
    }
  }
  return(hits)
}

modFindOverlap <- function(query_unlisted,subject_unlisted,subject_togroup,lenQuery, lenSubject, out){
  library(data.table)
  
  df <- cbind(start(subject_unlisted), end(subject_unlisted))
  options(digits=12)
  uniqueVector <- as.vector(matrix(df, ncol = 2) %*% matrix(rnorm(2),nrow = 2))
  posiciones_repetidas <- split(seq_along(uniqueVector), uniqueVector)[as.character(unique(uniqueVector))]
  new_subject_unlist <- df[!duplicated(uniqueVector),]
  rm(df)
  rm(uniqueVector)
  if (length(nrow(new_subject_unlist)) == 0){
    if (length(new_subject_unlist) == 2) {
      new_subject_unlist <- IRanges(start = new_subject_unlist[1], end = new_subject_unlist[2])
    }else{
      new_subject_unlist <- IRanges(start = NULL, end = NULL)
    } 
  }else{
    new_subject_unlist <- IRanges(start = new_subject_unlist[,1], end = new_subject_unlist[,2])
  }
  
  hits_listed_res <- list()
  positionsVector <- split(c(1:lenQuery), ceiling(seq_along(query_unlisted)/3000))
  for(val in c(1:length(positionsVector))){
    pos <- positionsVector[val]
    hits_listed <- as.list(findOverlaps(query_unlisted[unlist(pos)], new_subject_unlist, type = "any"))
    for (overlap in c(1:length(hits_listed))) {
      hits_listed[[overlap]] <- unique(subject_togroup[unlist(posiciones_repetidas[hits_listed[[overlap]]])])
    }
    if (length(hits_listed) == 0) {
      hits_listed_res <- hits_listed
      rm(hits_listed)
    }else{
      hits_listed_res <- c(hits_listed_res,hits_listed)
      rm(hits_listed)  
    }
    if(val%%10 == 0){
      gc()
    }
  }
  rm(posiciones_repetidas)
  gc()
  if (out != "list") {
    hits_listed_res <- Hits(rep(c(1:length(query_unlisted)),lengths(hits_listed_res)),unlist(hits_listed_res),lenQuery, lenSubject, sort.by.query = FALSE)
  }
  
  rm(new_subject_unlist)
  rm(subject_togroup)
  rm(query_unlisted)
  rm(subject_unlisted)
  
  gc()
  return(hits_listed_res)
}






splicesiteCounts <- function(x, frag_exonic, frag_intron, min_anchor,
                             option = c("junction", "exon"), include)
{
  
  option <- match.arg(option)
  
  N_L <- splicesiteOverlap(flank(x, -1, TRUE),
                           switch(option, junction = "R", exon = "L"),
                           frag_exonic, frag_intron, min_anchor, include)
  N_R <- splicesiteOverlap(flank(x, -1, FALSE),
                           switch(option, junction = "L", exon = "R"),
                           frag_exonic, frag_intron, min_anchor, include)
  N <- IntegerList(mapply(c, N_L, N_R, SIMPLIFY = FALSE))
  rm(frag_exonic)
  rm(frag_intron)
  rm(N_L)
  rm(N_R)
  gc(full=TRUE)
  return(N)
  
}

exonCoverage <- function(exons, exons_i_frag, frag_exonic)
{
  
  expanded_exon <- factor(togroup0(exons_i_frag), seq_along(exons))
  expanded_frag_exonic <- frag_exonic[unlist(exons_i_frag)]
  rm(exons_i_frag)
  rm(frag_exonic)
  
  expanded_exon <- expanded_exon[togroup0(expanded_frag_exonic)]
  expanded_frag_exonic <- unlist(expanded_frag_exonic)
  
  irl <- split(expanded_frag_exonic, expanded_exon)
  rm(expanded_exon)
  rm(expanded_frag_exonic)
  gc(full=TRUE)
  coverage <- coverage(irl, shift = -start(exons) + 1, width = width(exons))
  ## coverage() returns a SimpleRleList
  ## need RleList() to obtain a CompressedRleList
  coverage <- RleList(coverage)
  
  rm(exons)
  gc(full=TRUE)
  return(coverage)
  
}
