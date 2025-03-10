EventPointerStats_BAM <- function(PSI_boots, 
                                  Design, 
                                  Contrast, 
                                  Threshold = 0, 
                                  nbootstraps = 1000,
                                  cores=1,
                                  ram = 0.1,
                                  pathResult = "./"){

  if(is.null(Design)){
    stop("Design field is empty")
  }
  if(is.null(Contrast)){
    stop("Contrast field is empty")
  }
  pathResult <- paste0(pathResult, "/EventPointerStatsResult/")
  
  if (!file.exists(pathResult)) {
    dir.create(pathResult)
  } 
  
  result <- checkContrastDesignMatrices(Contrast, Design)
  if (result == TRUE){
    UseBootstrap <- T
    resBootstrap <- EventPointer_Bootstraps(PSI=PSI_boots, Design=Design,
                                            Contrast=Contrast,nBootstraps=nbootstraps,
                                            UsePseudoAligBootstrap =T,
                                            Threshold =Threshold,
                                            cores=cores, ram=ram)
    pathResultBootstrap <- paste0(pathResult, "bootstrapResult/")
    dir.create(pathResultBootstrap)
    for (coef in c(1:dim(resBootstrap$Pvalues)[2])){
      tableRes <- ResulTable(resBootstrap, coef = coef)
      write.csv(tableRes,file = paste0(pathResultBootstrap,"ResBootstrapContrast",coef,".csv"))
    }

  }
  
  
}

voomEventPointerBAM <- function(PSI_boots,Events,Design,Contrast){
  # Compute the abundance as the minimum value of the three paths
  # Using the minimum
  abundance <- unlist(sapply(Events,
                             function(y) sapply(y, function(x) min(rowMeans2(x$Counts)))))
  averageRef <- unlist(sapply(Events,
                              function(y) sapply(y, function(x) rowMeans2(x$Counts)[3])))
  averageP1 <- unlist(sapply(Events,
                             function(y) sapply(y, function(x) rowMeans2(x$Counts)[1])))
  averageP2 <- unlist(sapply(Events,
                             function(y) sapply(y, function(x) rowMeans2(x$Counts)[2])))

  # Remove some values
  dummy <- (rowSds(PSI_boots[,1,],useNames =T )<1e-6)
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


voom2 <- function (counts, design = NULL, lib.size = NULL, normalize.method = "none",
                   block = NULL, correlation = NULL, weights = NULL, span = 0.5,
                   plot = FALSE, save.plot = FALSE, keepMax=T)
{
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >0)
      design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size))
      lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts)))
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts)))
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  n <- nrow(counts)
  if (n < 2L)
    stop("Need at least two genes to fit a mean-variance trend")
  m <- min(counts)
  if (is.na(m))
    stop("NA counts not allowed")
  if (m < 0)
    stop("Negative counts not allowed")
  if (is.null(design)) {
    design <- matrix(1, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  if (is.null(lib.size))
    lib.size <- colSums(counts)
  y <- t(log2(t(counts + 0.5)/(lib.size + 1) * 1e+06))
  y <- limma:::normalizeBetweenArrays(y, method = normalize.method)
  y <- counts
  fit <- lmFit(y, design, block = block, correlation = correlation,
               weights = weights)
  if (is.null(fit$Amean))
    fit$Amean <- rowMeans(y, na.rm = TRUE)
  NWithReps <- sum(fit$df.residual > 0L)
  if (NWithReps < 2L) {
    if (NWithReps == 0L)
      warning("The experimental design has no replication. Setting weights to 1.")
    if (NWithReps == 1L)
      warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    out$design <- design
    if (is.null(out$targets))
      out$targets <- data.frame(lib.size = lib.size)
    else out$targets$lib.size <- lib.size
    return(new("EList", out))
  }
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e+06)
  sy <- sqrt(fit$sigma)
  allzero <- rowSums(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }
  l <- lowess(sx, sy, f = span)
  if(keepMax) {
    I <- which.max(l$y)
    ymax <- l$y[I]
    l$y[1:I] <- ymax
  }
  if (plot) {
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
         pch = 16, cex = 0.25)
    title("voom: Mean-variance trend")
    lines(l, col = "red")
  }
  f <- approxfun(l, rule = 2, ties = list("ordered", mean))
  if (fit$rank < ncol(design)) {
    j <- fit$pivot[1:fit$rank]
    fitted.values <- fit$coefficients[, j, drop = FALSE] %*%
      t(fit$design[, j, drop = FALSE])
  }
  else {
    fitted.values <- fit$coefficients %*% t(fit$design)
  }
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)
  w <- 1/f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)
  out$E <- y
  out$weights <- w
  out$design <- design
  if (is.null(out$targets))
    out$targets <- data.frame(lib.size = lib.size)
  else out$targets$lib.size <- lib.size
  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )",
                        ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }
  new("EList", out)
}
