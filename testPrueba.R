data(SG_RNASeq)
TxtPath<-tempdir()
PathSamplesAbundance <- system.file("extdata/bams", package = "EventPointer")
PathSamplesAbundance <- "D:/EventPointerBioconductor/EventPointer/inst/extdata/bams"

PathTranscriptomeGTF <- list.files(PathSamplesAbundance,"*.gtf",full.names = T)
PathTranscriptomeGTF <- "D:/EventPointer_3.0_replicate/simulation_data/annotation/refseq_hg19.formatted.gtf"
# TxDb <- GenomicFeatures:::makeTxDbFromGFF(file = PathTranscriptomeGTF,
#                                           format = "gtf", dataSource = "External Transcriptome")
PathSGResult <- "D:/EventPointerBioconductor/PruebaPaquete/"


BamInfo<-si
Samples<-BamInfo[,2]
library(EventPointer)
PathSamplesAbundance <- system.file("extdata/bams", package = "SGSeq")
PathTranscriptomeGTF<-paste(system.file("extdata",package="EventPointer"),"/FBXO31.gtf",sep="")
EventsDetection_BAM(PathSamplesAbundance, PathTranscriptomeGTF, 
                                      cores = 9, AnnEvents=T,
                                      min_junction_count = 5, max_complexity = 50,
                                      PathSGResult = PathSGResult)

load(file=paste0(PathSGResult,"PSI"))
Design <- cbind(rep(1,9),
                 rep(c(1,0,0),3),
                 rep(c(0,1,0),3))
#rownames(Dmatrix) <- nameBAM
Contrast <- cbind(
  c(0,1,0),
  c(0,0,1))

EventPointerStats(PSI_boots, Design, Contrast, cores=8, UseBootstrap=T,
                  Threshold = 0, ram = 0.1,
                  nbootstraps = 1000, pathResult = PathSGResult)
