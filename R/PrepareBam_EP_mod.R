
EventsDetection_BAM <- function(PathSamplesAbundance, 
                          PathTranscriptomeGTF = NULL, 
                          cores = 1, region=NULL,
                          min_junction_count = 2, max_complexity = 30,
                          lambda = NULL,nboot = 20,min_n_sample = NULL,
                          min_anchor = 6,
                          PathSGResult = ".") {
  
  
  cat("Getting BAM general information.../n")
  
  Bam_Info <- getBamInfo(PathSamplesAbundance,region = region, cores = cores)
  
  cat("Obtaining Reference Transcriptome...")

  
  TxDb <- GenomicFeatures::makeTxDbFromGFF(file = PathTranscriptomeGTF,
                                            format = "gtf", dataSource = "External Transcriptome")
  
  TxF_Ref <- convertToTxFeatures(TxDb)
  
  cat("Done")
  
  cat("\n Predicting Features from BAMs...")
  if (is.null(region)) {
    si <- seqinfo(BamFile(Bam_Info$file_bam[1]))
    sl <- rep(seqlevels(si),2)
    st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
    which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
    cat("\n Using this regions:\n")
    print(data.frame(which@seqnames,which@ranges))
  }else{
    si <- seqinfo(BamFile(Bam_Info$file_bam[1]))
    sl <- rep(seqlevels(si),2)
    st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
    which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
    which <- which[region]
    cat("\n Using this regions:\n")
    print(data.frame(which@seqnames,which@ranges))
    
  }
  
  if(is.null(min_n_sample)){
    min_n_sample<-min(c(length(Bam_Info$file_bam),3))
  }
  
  seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(which)
  
  TxF_RefLevels <- TxF_Ref[which(as.vector(seqnames(TxF_Ref)) %in% as.vector(seqnames(which)))]

  cat("\n Creating the splicing graph from the alignment files. This will take some time...")
  if (!is.null(which)) {
    TxF_mod <- predictTxFeatures(Bam_Info, cores = cores, which = which, 
                                 min_junction_count = min_junction_count,
                                 min_n_sample = min_n_sample, 
                                 max_complexity = max_complexity,
                                 min_anchor = min_anchor)  
  }else{
    TxF_mod <- predictTxFeatures(Bam_Info, cores = cores, 
                                 min_junction_count = min_junction_count,
                                 min_n_sample = min_n_sample,
                                 max_complexity = max_complexity,
                                 min_anchor = min_anchor)  
  }
  closeAllConnections()
  seqlevelsStyle(TxF_Ref) <- seqlevelsStyle(TxF_mod)
  
  features <- convertToSGFeatures(TxF_mod)
  features <- annotate(features, TxF_RefLevels)
  valid_rows <- which(lengths(features@geneName)>0)
  features <- features[valid_rows]
  cat("\n Assigning read counts to the paths in the splicing graph...")
  SgFC <- getSGFeatureCounts(Bam_Info, min_anchor = min_anchor,
                                                features, cores = cores)
  closeAllConnections()
  cat("\n Detect event from splicing graph...")
  SgFC <- annotate(SgFC, TxF_RefLevels)
  
  save(SgFC,file=paste0(PathSGResult,"/SgFC.RData"))
  EventsDetection_pred<-EventDetection(SgFC, cores=cores)
  
  cat("\n Checking that detected events are previously annotated...")

  save(EventsDetection_pred, file=paste0(PathSGResult,"/EventsDetection_EPBAM.RData"))
  
  PSI_boots <- getPSI_RNASeq_boot(EventsDetection_pred, lambda = lambda, cores, nboot)
  save(PSI_boots, file=paste0(PathSGResult,"/PSI_boots.RData"))
  
  event_list <- list()
  index <- 1
  
  for (geneEvent in EventsDetection_pred) {
    for (event in geneEvent) {
      event_list[[index]] <- event$Info
      index <- index + 1
    }
  }
  
  totalEventTable <- do.call(rbind, event_list)
  
  write.csv(totalEventTable, file = paste0(PathSGResult, "/TotalEventsFound.csv"), row.names = FALSE)
  
  cat("\n DONE!")
  closeAllConnections()

}

AnnEventsFunc <- function(EventsDetection_pred, EventsDetection_ann, cores){
  library(doParallel)
  registerDoParallel(cores=cores)
  listEventsPred <- foreach(event=unlist(EventsDetection_pred, recursive = F), .packages = 'GenomicRanges') %dopar% {
    P1 <- GRanges(event$P1)
    P2 <- GRanges(event$P2)
    Ref <- GRanges(event$Ref)
    GRangesList(c(P1,P2,Ref))[[1]]
  }
  listEventsPred <- GRangesList(listEventsPred)
  # listEventsAnn <- GRangesList()
  listEventsAnn <- foreach(event=unlist(EventsDetection_ann, recursive = F), .packages = 'GenomicRanges') %dopar% {
    P1 <- GRanges(event$P1)
    P2 <- GRanges(event$P2)
    Ref <- GRanges(event$Ref)
    GRangesList(c(P1,P2,Ref))[[1]]
  }
  stopImplicitCluster()
  listEventsAnn <- GRangesList(listEventsAnn)
  fso <- GenomicAlignments:::findSpliceOverlaps(listEventsPred,listEventsAnn, type="equal")
  fso <- fso[which(fso@elementMetadata$compatible)]
  fso <- fso[which(lengths(listEventsAnn[fso@to])==lengths(listEventsPred[fso@from]))]
  percentIdentity <- 1-sum(abs(width(listEventsAnn[fso@to])-width(listEventsPred[fso@from])))/max(sum(width(listEventsAnn[fso@to])), sum(width(listEventsPred[fso@from])))
  fso <- fso[which(percentIdentity >=.98)]
  
  countEvent <- 1
  for (gene in c(1:length(EventsDetection_pred))) {
    for(posEvent in c(1:length(EventsDetection_pred[[gene]]))){
      event <- EventsDetection_pred[[gene]][[posEvent]]
      if(countEvent %in% fso@from){
        event$Info$Ann <- T
      } else{
        event$Info$Ann <- F
      }
      EventsDetection_pred[[gene]][[posEvent]] <- event
      countEvent <- countEvent+1 
    }
          
  }
  return(EventsDetection_pred)
}

getBamInfo <- function(PathSamplesAbundance, region, cores = 1)
{
  sample_name <- dir(PathSamplesAbundance,pattern = "*.bam$")
  if(identical(sample_name,character(0))){
    file_bam <- c()
    for(folder in list.dirs(PathSamplesAbundance)){
      files_bam_folder <- list.files(folder,pattern = "*.bam$")
      sample_name <- c(sample_name,files_bam_folder)
      file_bam <- c(file_bam, rep(folder, length(files_bam_folder)))
      
    }
    sample_info <- data.frame(sample_name = sample_name,
                              file_bam = file_bam, stringsAsFactors = FALSE)
  }else{
    sample_info <- data.frame(sample_name = sample_name,
                              file_bam = rep(PathSamplesAbundance, length(sample_name)), stringsAsFactors = FALSE)
  }
  
  SGSeq:::checkSampleInfo(sample_info, FALSE)
  listSamples <- split(sample_info, c(1:nrow(sample_info)))
  clF <- makePSOCKcluster(cores, outfile = "bamInfoSamples.txt")
  list_bamInfo <- clusterApplyLB(cl = clF,
                                 x = listSamples,
                                 fun = getBamInfoPerSample,
                                 region=region)
  on.exit(stopCluster(cl = clF))
  SGSeq:::checkApplyResultsForErrors(list_bamInfo,
                             "getBamInfoPerSample",
                             sample_info$sample_name,
                             "character")
  
  bamInfo <- do.call(SGSeq:::rbindDfsWithoutRowNames, list_bamInfo)
  
  SGSeq:::checkBamInfo(bamInfo)
  
  cols <- c("paired_end", "read_length", "frag_length", "lib_size")
  cols <- cols[cols %in% names(bamInfo)]
  
  for (col in cols) {
    sample_info[[col]] <- bamInfo[[col]]
    
  }
  sample_info$file_bam <- paste0(sample_info$file_bam,"/",sample_info$sample_name)
  return(sample_info)
  
}


getBamInfoPerSample <- function(sample_info, region)
{
  file_bam <- paste0(sample_info$file_bam,"/",sample_info$sample_name)
  sample_name <- sample_info$sample_name
  if (is(file_bam, "BamFile")) {
    file_tmp <- file_bam
    
  } else {
    file_tmp <- BamFile(file_bam)
    
  }
  
  si <- seqinfo(file_tmp)
  sl <- rep(seqlevels(si),2)
  st <- c(rep(c("+"), length(si)),rep(c("-"), length(si)))
  which <- GRanges(sl, IRanges(1, seqlengths(si)[sl]), st)
  
  flag <-
    scanBamFlag(
      isUnmappedQuery = FALSE,
      isSecondaryAlignment = FALSE,
      hasUnmappedMate = FALSE,
      isDuplicate = FALSE
      
    )
  what <- c("qname", "flag", "qwidth", "isize")
  param <-
    ScanBamParam(
      flag = flag,
      what = what,
      which = which[region[1]],
      tag = "XS"
    )
  
  bam <- scanBam(file = file_tmp, param = param)[[1]]
  gc()
  XS <- !is.null(bam$tag$XS)
  paired_end <- any(bamFlagTest(bam$flag, "isPaired"))
  read_length <- median(bam$qwidth, na.rm = TRUE)
  x <- data.frame(
    XS = XS,
    paired_end = paired_end,
    read_length = read_length,
    stringsAsFactors = FALSE
  )
  if (paired_end) {
    isize <- bam$isize
    frag_length <- median(isize[which(isize > 0)], na.rm = TRUE)
    infoBam <- idxstatsBam(file = file_bam, index = file_bam)
    x$lib_size <- sum(infoBam$mapped)/2
    
  } else {
    infoBam <- idxstatsBam(file = file_bam, index = file_bam)
    desv <- (infoBam$mapped[1]/length(unique(bam$qname))-1)
    x$lib_size <- sum(infoBam$mapped) - (sum(infoBam$mapped) * desv)
    frag_length <- NA_real_
    
  }
  rm(bam)
  gc()
  x$frag_length <- frag_length
  
  
  message(sample_name)
  
  return(x)
  
}
