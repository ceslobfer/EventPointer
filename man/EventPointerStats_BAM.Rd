% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Prepare_EP_mod.R
\name{EventPointerStats_BAM}
\alias{EventPointerStats_BAM}
\title{Differential comparison between sample groups for the calculation of delta PSI}
\usage{
EventPointerStats_BAM(PSI_boots, 
                      Design, 
                      Contrast, 
                      Threshold = 0,
                      nbootstraps = 1000,
                      cores=1, 
                      ram = 0.1, 
                      pathResult = PathSGResult)
}

\description{
This function performs the differential comparison between sample groups for the calculation of delta PSI from \code{EventDetection_BAM} function.
}

\arguments{

\item{PSI_boots}{PSI_boots.RData obtained in \code{EventsDetection_BAM} function.}

\item{Design}{A matrix defining the linear model. Each row corresponds to an array, and each column corresponds to a coefficient (such as the baseline and treatment effects).}

\item{Contrast}{A numeric matrix with contrasts to be tested. Rows correspond to coefficients in the design matrix, and columns correspond to contrasts.}

\item{Threshold}{It assigns a threshold to compute the pvalues. Default value is 0.}

\item{nbootstraps}{How many layers, Bootstraps or samplings are goingto be used. Caution, high numbers increase computational time.}

\item{cores}{Number of cores used for parallel processing}

\item{pathResult}{Folder where save the function result.}
}
\value{
Folder containing .csv files with the results associated with the differential \eqn{\Psi} analysis of the indicated contrasts.It presents the \eqn{\Delta \Psi} associated with each event of that contrast and its corresponding significance parameters. Note that a table will be obtained for each contrast indicated in the contrast matrix.

}

\examples{
\dontrun{
  # Run EventsDetection_BAM function
   data(PSI_boots)
   
   Design <- cbind(rep(1,9),rep(c(1,0,0),3),rep(c(0,1,0),3))
   Contrast <- cbind(c(0,1,0),c(0,0,1))
   
   EventPointerStats_BAM(PSI_boots, 
                      Design, 
                      Contrast, 
                      Threshold = 0,
                      nbootstraps = 1000,
                      cores=1, 
                      ram = 0.1, 
                      pathResult = PathSGResult)
}
}