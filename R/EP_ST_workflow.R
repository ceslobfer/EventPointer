
#General first step for EventPointer ST workflow
EventsDetection_ST<- function(PathSamplesAbundance=NULL,
                              typeAbundance = "kallisto", 
                              PathTranscriptomeGTF = NULL,
                              EventsTranscriptome=NULL,
                              Bootstrap=T,
                              Filter=F,
                              Qn = 0.25,
                              cores=1,
                              PathEventsGTFResults="."){

  #Event detection from transcriptome gtf
  if(!is.null(PathTranscriptomeGTF) & is.null(EventsTranscriptome)){ 
    EventsTranscriptome <- EventDetection_transcriptome(inputFile = PathTranscriptomeGTF,
                                                        Pathtxt=PathEventsGTFResults,
                                                        cores=cores)
    save(EventsTranscriptome,file=paste0(PathEventsGTFResults,"/EventsTranscriptome.RData"))
  }
  
  
  
  if(!is.null(EventsTranscriptome) & !is.null(PathSamplesAbundance)){ 
    data_exp <- getbootstrapdata(PathSamples = PathSamplesAbundance,type = typeAbundance)
    
    rownames(data_exp[[1]]) <- sapply(strsplit(rownames(data_exp[[1]]),"\\|"),function(X) return(X[1]))
    
    PSI <- GetPSI_FromTranRef(PathsxTranscript = EventsTranscriptome,
                              Samples = data_exp,
                              Bootstrap = Bootstrap,
                              Filter = Filter,
                              Qn = Qn)
  }else{
    stop("Events not provided")
  }
  save(PSI,file=paste0(PathEventsGTFResults,"/PSI_ST.RData"))
  return(PSI)
  
  
}

EventPointerStats_ST <- function(PSI, 
                                 Design, 
                                 Contrast, 
                                 BootstrapStats = T,
                                 nbootstraps= 10000, 
                                 UsePseudoAligBootstrap = T,
                                 Threshold = 0,
                                 cores=1, 
                                 ram=4,
                                 pathResult="./"){
  
  if(is.null(Design)){
    stop("Design field is empty")
  }
  if(is.null(Contrast)){
    stop("Contrast field is empty")
  }
  pathResult <- paste0(pathResult, "/EventPointerStatsSTResult/")
  checkMatrices <- checkContrastDesignMatrices(Contrast, Design)
  if(checkMatrices){
    if (!file.exists(pathResult)) {
      dir.create(pathResult)
    } 
    
    if(BootstrapStats){
      resBootstrap <- EventPointer_Bootstraps(PSI$PSI,
                                              Design,
                                              Contrast,
                                              cores = cores,
                                              ram = ram,
                                              nbootstraps = nbootstraps,
                                              UsePseudoAligBootstrap = UsePseudoAligBootstrap,
                                              Threshold = Threshold)
      pathResultBootstrap <- paste0(pathResult, "bootstrapResult/")
      dir.create(pathResultBootstrap)
      for (coef in c(1:dim(resBootstrap$Pvalues)[2])){
        tableRes <- ResulTable(resBootstrap, coef = coef)
        write.csv(tableRes,file = paste0(pathResultBootstrap,"ResBootstrapContrast",coef,".csv"))
      }
    }
    
  }
  
}


voomEventPointerST <- function(PSI,Design,Contrast){
  # Compute the abundance as the minimum value of the three paths
  # Using the minimum
  PSI_boots <- PSI$PSI
  Events <- PSI$ExpEvs
  abundance <- unlist(lapply(Events,
                             function(y) min(colMeans2(y))))
  
  averageP1 <- unlist(lapply(Events,
                             function(y) colMeans2(y)[1]))
  averageP2 <- unlist(lapply(Events,
                             function(y) colMeans2(y)[2]))
  averageRef <- unlist(lapply(Events,
                              function(y) colMeans2(y)[3]))
  
  # Remove some values
  dummy <- (rowSds(PSI_boots[,1,],useNames =T)<1e-6)
  dummy[is.na(dummy)] <- TRUE
  Quitar <- dummy
  
  PSI <- PSI_boots[!Quitar,1,]
  abundance <- abundance[!Quitar]
  averageP1 <- averageP1[!Quitar]
  averageP2 <- averageP2[!Quitar]
  averageRef <- averageRef[!Quitar]
  betaest <- solve(base::crossprod(Design),t(Design)%*%t(PSI))
  error <- t(PSI) - Design %*% betaest
  errorSE <- sqrt(colSums(error^2)/-diff(dim(Design)))
  
  # If the value of PSI is close to 1 or close to zero, the variance is small
  p <- rowMeans(PSI)
  varEst <- p*(1-p)
  # Let's combine both sources of information
  modelo <- lm((errorSE)^.5 ~ I(log2(abundance)) +
                 I(log2(averageRef)) + I(log2(averageP1)) + I(log2(averageP2)))
  
  overallVar <- predict(modelo)
  
  overallVar <- 1 + max(overallVar) - overallVar
  
  logabd <- overallVar
  PSI_changed <- (5+PSI) * logabd/rowMeans(5+PSI)
  
  # We will try to apply voom to these data.
  # Remove NAs for the time being
  
  modelo <- voom2(PSI_changed, Design, plot=F, span=.1)
  fit <- suppressWarnings(lmFit(modelo, Design, method="robust"))
  fit <- contrasts.fit(fit, Contrast)
  fit <- suppressWarnings(eBayes(fit, robust=TRUE,proportion = .1))
  deltaPSI <- t(Contrast) %*% (solve(base::crossprod(Design), t(Design))) %*% t(PSI)
  colnames(deltaPSI) <- rownames(PSI_changed)
  res <- list()
  for(coef in c(1:ncol(fit$coefficients))){
    result <- topTable(fit, coef = coef, number = Inf)
    result <- data.frame(deltaPSI = deltaPSI[coef,rownames(result)], result)
    res[[coef]] <- result
  }
  return(res)
}
