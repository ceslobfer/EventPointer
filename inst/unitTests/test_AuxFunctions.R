test_AuxFunctions <- function() {
  
  checkIdentical(estimateAbsoluteConc(50,50,100,1)$T1est,100)
  checkIdentical(estimateAbsoluteConc(25,50,100,10)$T2est,100)
  
  Mat<-matrix(c(0,-1,1,-1,-1,1,1,-1,1,1),ncol=2)
  checkEqualsNumeric(getRandomFlow(Mat)[2,1],-1,tolerance=1e-8)
  
  checkEqualsNumeric(IHsummarization(0.5,2,0.001,10)$Pvalues,0.5632502,tolerance=0.001)
  
  checkEqualsNumeric(sum(pdist2(diag(10),diag(10))),180,tolerance=0.001)
  
}
