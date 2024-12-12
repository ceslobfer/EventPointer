library(EventPointer)
data(SG_RNASeq)
TxtPath<-tempdir()
PathSamplesAbundance <- system.file("extdata/bams", package = "EventPointer")
PathSamplesAbundance <- "D:/EventPointerBioconductor/EventPointer/inst/extdata/bams"

PathTranscriptomeGTF <- list.files(PathSamplesAbundance,"*.gtf",full.names = T)
PathTranscriptomeGTF <- "D:/EventPointer_3.0_replicate/simulation_data/annotation/refseq_hg19.formatted.gtf"

PathSGResult <- "D:/EventPointerToolkit/PruebaPaquete/"


BamInfo<-si
Samples<-BamInfo[,2]
library(EventPointer)
PathSamplesAbundance <- system.file("extdata/bams", package = "SGSeq")
PathTranscriptomeGTF<-paste(system.file("extdata",package="EventPointer"),"/FBXO31.gtf",sep="")
EventsDetection_BAM(PathSamplesAbundance, PathTranscriptomeGTF, 
                                      cores = 9,
                                      min_junction_count = 5, max_complexity = 50,
                                      PathSGResult = PathSGResult)

load(file=paste0(PathSGResult,"PSI"))

Design <- cbind(rep(1,9),
                 rep(c(1,0,0),3),
                 rep(c(0,1,0),3))
Contrast <- cbind(
  c(0,1,0),
  c(0,0,1))

EventPointerStats_BAM(PSI_boots, Design, Contrast, cores=8, pathResult = PathSGResult)

EventDetection()
SG_RNASeq <- load("D:/EventPointerBioconductor/PruebaPaquete/SgFC.RData")
SgFC
load("D:/EventPointerBioconductor/PruebaPaquete/PSI_boots.RData")
EventsTxt <- read.csv("D:/EventPointerBioconductor/PruebaPaquete/TotalEventsFound.csv",row.names = 1)
EventsTxt<-"D:/EventPointerBioconductor/PruebaPaquete/TotalEventsFound.csv"
PathGTF <- "D:/EventPointer/vignettes/"
EventPointerBAM_IGV(SG_RNASeq, EventsTxt, PathGTF)


EventsDetection_ST(PathSamplesAbundance = NULL, PathTranscriptomeGTF = NULL, EventsTranscriptome = NULL,
                   PathEventsGTFResults=".",
                   cores=1, typeAbundance = "kallisto", Bootstrap=T, Filter=F,
                   Qn = 0.25)

PathFiles<-system.file("extdata",package="EventPointer")
PathTranscriptomeGTF <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep="")
PathSamplesAbundance <- paste0(PathFiles,"/output")
PathSamplesAbundance <- dir(PathSamplesAbundance,full.names = TRUE)
Pathtxt <- "D:/EventPointerBioconductor/PruebaPaquete/"
EventsPSI <- EventsDetection_ST(PathSamplesAbundance,PathTranscriptomeGTF = PathTranscriptomeGTF,
                                PathEventsGTFResults=Pathtxt,
                                cores=1, typeAbundance = "kallisto", Bootstrap=T, Filter=F,
                                Qn = 0.25)
Design <- cbind(1,rep(c(0,1),each=2))
Contrast <- matrix(c(0,1),nrow=2)
load("D:/EventPointerBioconductor/PruebaPaquete/EventsTranscriptome.RData")
EventsTranscriptome

EventPointerStats_ST(EventsPSI, Design, Contrast, cores=1, ram=4, 
                                 BootstrapStats = T,nbootstraps= 10000, 
                                 UsePseudoAligBootstrap = T,Threshold = 0,
                                 pathResult="D:/EventPointerBioconductor/PruebaPaquete/")

