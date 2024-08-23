## ----LoadFunctions, echo=FALSE, message=FALSE, warning=FALSE, results='hide'----
library(knitr)
opts_chunk$set(error = FALSE)
library(EventPointer)
library(dplyr)
library(kableExtra)

## ----style, echo = FALSE, results = 'asis'------------------------------------
##BiocStyle::markdown()

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(BiocManager)
#  
#  if (!requireNamespace("BiocManager", quietly=TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("EventPointer")

## ----CDFGTF, eval=FALSE, warning=FALSE, collapse=TRUE-------------------------
#  
#  # Set input variables
#     PathFiles<-system.file("extdata",package="EventPointer")
#     DONSON_GTF<-paste(PathFiles,"/DONSON.gtf",sep="")
#     PSRProbes<-paste(PathFiles,"/PSR_Probes.txt",sep="")
#     JunctionProbes<-paste(PathFiles,"/Junction_Probes.txt",sep="")
#     Directory<-tempdir()
#     array<-"HTA-2_0"
#  
#  # Run the function
#  
#     CDFfromGTF(input="AffyGTF",inputFile=DONSON_GTF,
#                PSR=PSRProbes,Junc=JunctionProbes,
#                PathCDF=Directory,microarray=array)
#  

## ----aroma, eval=FALSE--------------------------------------------------------
#  
#  # Simple example of Aroma.Affymetrix Preprocessing Pipeline
#  
#  verbose <- Arguments$getVerbose(-8);
#  timestampOn(verbose);
#  projectName <- "Experiment"
#  cdfGFile <- "EP_HTA-2_0,r"
#  cdfG <- AffymetrixCdfFile$byChipType(cdfGFile)
#  cs <- AffymetrixCelSet$byName(projectName, cdf=cdfG)
#  bc <- NormExpBackgroundCorrection(cs, method="mle", tag=c("*","r11"));
#  csBC <- process(bc,verbose=verbose,ram=20);
#  qn <- QuantileNormalization(csBC, typesToUpdate="pm");
#  csN <- process(qn,verbose=verbose,ram=20);
#  plmEx <- ExonRmaPlm(csN, mergeGroups=FALSE)
#  fit(plmEx, verbose=verbose)
#  cesEx <- getChipEffectSet(plmEx)
#  ExFit <- extractDataFrame(cesEx, addNames = TRUE)

## ----EP_arrays, eval=TRUE-----------------------------------------------------

   data(ArraysData)

   Dmatrix<-matrix(c(1,1,1,1,0,0,1,1),nrow=4,ncol=2,byrow=FALSE)
   Cmatrix<-t(t(c(0,1)))
   EventsFound<-paste(system.file("extdata",package="EventPointer"),"/EventsFound.txt",sep="")
   
   Events<-EventPointer(Design=Dmatrix,
                      Contrast=Cmatrix,
                      ExFit=ArraysData,
                      Eventstxt=EventsFound,
                      Filter=FALSE,
                      Qn=0.25,
                      Statistic="LogFC",
                      PSI=TRUE)

## ----EP_Arrays_Res_Table, echo=FALSE------------------------------------------
kable(Events[1:5,],digits=5,row.names=TRUE,align="c",caption = "Table 1: EventPointer Arrays results")

## ----Arrays_IGV, eval=FALSE, collapse=TRUE------------------------------------
#  
#  # Set Input Variables
#  
#     DONSON_GTF<-paste(PathFiles,"/DONSON.gtf",sep="")
#     PSRProbes<-paste(PathFiles,"/PSR_Probes.txt",sep="")
#     JunctionProbes<-paste(PathFiles,"/Junction_Probes.txt",sep="")
#     Directory<-tempdir()
#     EventsFound<-paste(system.file("extdata",package="EventPointer"),"/EventsFound.txt",sep="")
#     array<-"HTA-2_0"
#  
#  
#  # Generate Visualization files
#  
#   EventPointer_IGV(Events[1,,drop=FALSE],"AffyGTF",DONSON_GTF,PSRProbes,JunctionProbes,Directory,EventsFound,array)
#  
#  

## ----PrepareBam, eval=FALSE, collapse=TRUE------------------------------------
#  # Obtain the samples and directory for .bam files
#  
#  # the object si contains example sample information from the SGSeq R package
#  # use ?si to see the corresponding documentation
#  
#     BamInfo<-si
#     Samples<-BamInfo[,2]
#     PathToSamples <- system.file("extdata/bams", package = "SGSeq")
#     PathToGTF<-paste(system.file("extdata",package="EventPointer"),"/FBXO31.gtf",sep="")
#  
#    # Run PrepareBam function
#     SG_RNASeq<-PrepareBam_EP(Samples=Samples,
#                              SamplePath=PathToSamples,
#                              Ref_Transc="GTF",
#                              fileTransc=PathToGTF,
#                              cores=1)

## ----EventDetection, eval=TRUE------------------------------------------------
  # Run EventDetection function
   data(SG_RNASeq)
   TxtPath<-tempdir()
   AllEvents_RNASeq<-EventDetection(SG_RNASeq,cores=1,Path=TxtPath)
   

## ----ListofLists, eval=FALSE--------------------------------------------------
#  Events[[i]][[j]]

## ----EP_RNASeq, eval=TRUE-----------------------------------------------------
   Dmatrix<-matrix(c(1,1,1,1,1,1,1,1,0,0,0,0,1,1,1,1),ncol=2,byrow=FALSE)
   Cmatrix<-t(t(c(0,1)))
   Events <- EventPointer_RNASeq(AllEvents_RNASeq,Dmatrix,Cmatrix,Statistic="LogFC",PSI=TRUE)

## ----EP_RNASeq_Res_Table, echo=FALSE------------------------------------------
kable(Events[1:5,],digits=5,row.names=TRUE,align="c",caption = "Table 2: EventPointer RNASeq results")

## ----RNAS_IGV, eval=TRUE, collapse=TRUE---------------------------------------

   # IGV Visualization
   EventsTxt<-paste(system.file("extdata",package="EventPointer"),"/EventsFound_RNASeq.txt",sep="")
   PathGTF<-tempdir()
   EventPointer_RNASeq_IGV(Events,SG_RNASeq,EventsTxt,PathGTF)

## ----GTFfGTF, eval=TRUE, warning=FALSE, collapse=TRUE-------------------------


# Set input variables
PathFiles<-system.file("extdata",package="EventPointer")
inputFile <- paste(PathFiles,"/gencode.v24.ann_2genes.gtf",sep="")
Transcriptome <- "Gencode24_2genes"
Pathtxt <- tempdir()

## ----GTFfGTF_2, eval=FALSE, warning=FALSE, collapse=TRUE----------------------
#  # Run the function
#  inputFile <- "/home/external_HDD/HDD_1/EventPointer/inst/extdata/gencode.v24.ann_2genes_prueba.gtf"
#  EventXtrans <- EventDetection_transcriptome(inputFile = inputFile,
#                                               Transcriptome = Transcriptome,
#                                               Pathtxt=Pathtxt,
#                                               cores=1)
#  

## ----GTFfGTFnames, eval=FALSE, warning=FALSE, collapse=TRUE-------------------
#  
#  names(EventXtrans)
#  

## ----getbootstrpdata, eval=TRUE, warning=FALSE, collapse=TRUE-----------------

PathSamples<-system.file("extdata",package="EventPointer")
PathSamples <- paste0(PathSamples,"/output")
PathSamples <- dir(PathSamples,full.names = TRUE)

data_exp <- getbootstrapdata(PathSamples = PathSamples,type = "kallisto")


## ----getdata_nb, eval=TRUE, warning=FALSE, collapse=TRUE----------------------

# this code chunk is an example of how to load the data from kallisto output.
# the expression of the isoforms are counts

PathFiles <- system.file("extdata",package="EventPointer")
filesnames <- dir(paste0(PathFiles,"/output"))
PathFiles <- dir(paste0(PathFiles,"/output"),full.names = TRUE)
dirtoload <- paste0(PathFiles,"/","abundance.tsv")
RNASeq <- read.delim(dirtoload[1],sep = "\t", colClasses = c(NA,"NULL","NULL",NA,"NULL"))
for (n in 2:length(dirtoload)){
  RNASeq[,n+1] <- read.delim(dirtoload[n],sep = '\t', colClasses = c('NULL','NULL','NULL',NA,'NULL'))
}
rownames(RNASeq)<-RNASeq[,1]
RNASeq<-RNASeq[,-1]
colnames(RNASeq) <- filesnames



## ----check_ann_nb, eval=TRUE, warning=FALSE, collapse=TRUE--------------------
rownames(RNASeq)[1:5]
EventXtrans$transcritnames[1:5]

## ----PSI_Statistic, eval=TRUE, warning=FALSE, collapse=TRUE-------------------

#change rownames of RNASeq variable
rownames(RNASeq) <- sapply(strsplit(rownames(RNASeq),"\\|"),function(X) return(X[1]))
RNASeq<-as.matrix(RNASeq) #must be a matrix variable


PSIss_nb <- GetPSI_FromTranRef(PathsxTranscript = EventXtrans,
                            Samples = RNASeq,
                            Bootstrap = FALSE,
                            Filter = FALSE)

PSI <- PSIss_nb$PSI
Expression <- PSIss_nb$ExpEvs


## ----check_ann, eval=TRUE, warning=FALSE, collapse=TRUE-----------------------
rownames(data_exp[[1]])[1:5]
EventXtrans$transcritnames[1:5]


## ----psi_with_bootstrap, eval=TRUE, warning=FALSE, collapse=TRUE--------------

#change rownames of the first element of teh list data_exp
rownames(data_exp[[1]]) <- sapply(strsplit(rownames(data_exp[[1]]),"\\|"),function(X) return(X[1]))

PSIss <- GetPSI_FromTranRef(PathsxTranscript = EventXtrans,
                               Samples = data_exp,
                               Bootstrap = TRUE,
                               Filter = FALSE)

PSI <- PSIss$PSI
Expression <- PSIss$ExpEvs

## ----PSI_Statistic2, eval=TRUE, warning=FALSE, collapse=TRUE------------------

# Design and contrast matrix:

Design <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
Contrast <- matrix(c(0,1),nrow=2)

# Statistical analysis:

Fit <- EventPointer_RNASeq_TranRef(Count_Matrix = Expression,Statistic = "LogFC",Design = Design, Contrast = Contrast)


## ----EP_TranRef_Res_Table, echo=FALSE-----------------------------------------
kable(Fit,digits=5,row.names=FALSE,align="c",caption = "Table 3: PSI_Statistic results")

## ----ep_bootstrp_statistic, eval=TRUE, warning=FALSE, collapse=TRUE-----------

Dmatrix <- cbind(1,rep(c(0,1),each=2))
Cmatrix <- matrix(c(0,1),nrow=2)

Fit <- EventPointer_Bootstraps(PSI = PSI,
                               Design = Dmatrix,
                               Contrast = Cmatrix,
                               cores = 1,
                               ram = 1,
                               nBootstraps = 10,
                               UsePseudoAligBootstrap = TRUE)


## ----ResultTable, eval=TRUE, warning=FALSE, collapse=TRUE---------------------

ResulTable(EP_Result = Fit,coef = 1,number = 5)


## ----ep_tranref_igv, eval=TRUE, warning=FALSE, collapse=TRUE------------------

SG_List <- EventXtrans$SG_List
PathEventsTxt<-system.file('extdata',package='EventPointer')
PathEventsTxt <- paste0(PathEventsTxt,"/EventsFound_Gencode24_2genes.txt")
PathGTF <- tempdir()

EventPointer_RNASeq_TranRef_IGV(SG_List = SG_List,pathtoeventstable = PathEventsTxt,PathGTF = PathGTF)


## ----biomart, eval=FALSE, warning=FALSE, collapse=TRUE------------------------
#  
#  library(biomaRt)
#  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "mar2016.archive.ensembl.org")
#  mart<-useDataset("hsapiens_gene_ensembl",mart)
#  
#  mistranscritos <- EventXtrans$transcritnames
#  head(mistranscritos)
#  
#  ## we need to remove the ".x":
#  mistranscritos <- gsub("\\..*","",mistranscritos)
#  
#  Dominios <- getBM(attributes = c("ensembl_transcript_id","interpro","interpro_description"),
#                    filters = "ensembl_transcript_id",
#                    values = mistranscritos,
#                    mart=mart)
#  
#  #we build the isoform x protein domain matrix
#  library(Matrix)
#  ii <- match(Dominios$ensembl_transcript_id,mistranscritos)
#  misDominios <- unique(Dominios[,2])
#  jj <- match(Dominios[,2],misDominios)
#  
#  TxD <- sparseMatrix(i=ii,j = jj,dims = c(length(mistranscritos),length(misDominios)),x = 1)
#  rownames(TxD) <- mistranscritos
#  colnames(TxD) <- misDominios
#  

## ----pd_enrichmetn, eval=FALSE, warning=FALSE, collapse=TRUE------------------
#  
#  data("TxD")
#  
#  #check same annotation for transcripts:
#  EventXtrans$transcritnames[1]
#  # [1] "ENST00000611540.4"
#  rownames(TxD)[1]
#  # [1] "ENST00000611540"
#  
#  
#  ## as is not the same, we change EventXtrans$transcritnames annotation
#  transcriptnames <- EventXtrans$transcritnames
#  transcriptnames <- gsub("\\..*","",transcriptnames)
#  EventXtrans$transcritnames <- transcriptnames
#  
#  Result_PDEA <- Protein_Domain_Enrichment(PathsxTranscript = EventXtrans,
#                                           TxD = TxD,
#                                           Diff_PSI = Fit$deltaPSI)
#  
#  

## ----pathprimer3, eval=FALSE, warning=FALSE, collapse=TRUE--------------------
#  
#  Primer3Path <- Sys.which("primer3_core")
#  

## ----FindPrimers, eval=FALSE, warning=FALSE, collapse=TRUE--------------------
#  
#  data("EventXtrans")
#  #From the output of EventsGTFfromTranscriptomeGTF we take the splicing graph information
#  SG_list <- EventXtrans$SG_List
#  #SG_list contains the information of the splicing graphs for each gene
#  #Let's supone we want to design primers for the event 1 of the gene ENSG00000254709.7
#  
#  #We take the splicing graph information of the required gene
#  SG <- SG_list$ENSG00000254709.7
#  
#  #We point the event number
#  EventNum <- 1
#  
#  #Define rest of variables:
#  Primer3Path <- Sys.which("primer3_core")
#  Dir <- "C:\\PROGRA~2\\primer3\\"
#  
#  
#  MyPrimers_taqman <- FindPrimers(SG = SG,
#                                  EventNum = EventNum,
#                                  Primer3Path = Primer3Path,
#                                  Dir = Dir,
#                                  mygenomesequence = BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
#                                  taqman = 1,
#                                  nProbes=1,
#                                  nPrimerstwo=4,
#                                  ncommonForward=4,
#                                  ncommonReverse=4,
#                                  nExons=10,
#                                  nPrimers =5,
#                                  maxLength = 1200)
#  

## ----EP_DesiCP_Res_Table, eval=TRUE, warning=FALSE, collapse=TRUE,echo=FALSE----
data("MyPrimers")
kable(MyPrimers[1:5,],digits=5,row.names=FALSE,align="c",caption = "Table 4: Data.frame output of FindPrimers for conventional PCR") %>%
  kable_styling() %>%
  scroll_box(width ="660px")


## ----EP_DesiTP_Res_Table, eval=TRUE, warning=FALSE, collapse=TRUE,echo=FALSE----
data("MyPrimers_taqman")
kable(MyPrimers_taqman[1:5,],digits=5,row.names=FALSE,align="c",caption = "Table 5: Data.frame output of FindPrimers for conventional PCR") %>%
  kable_styling() %>%
  scroll_box(width ="660px")


## ----CDFGTF_MP, eval=FALSE, warning=FALSE, collapse=TRUE----------------------
#  
#  # Set input variables
#     PathFiles<-system.file("extdata",package="EventPointer")
#     DONSON_GTF<-paste(PathFiles,"/DONSON.gtf",sep="")
#     PSRProbes<-paste(PathFiles,"/PSR_Probes.txt",sep="")
#     JunctionProbes<-paste(PathFiles,"/Junction_Probes.txt",sep="")
#     Directory<-tempdir()
#     array<-"HTA-2_0"
#  
#  # Run the function
#  
#     CDFfromGTF_Multipath(input="AffyGTF",inputFile=DONSON_GTF,
#                PSR=PSRProbes,Junc=JunctionProbes,
#                PathCDF=Directory,microarray=array,paths=3)
#  

## ----EventDetection_MP, eval=FALSE--------------------------------------------
#    # Run EventDetection function
#     data(SG_RNASeq)
#     TxtPath<-tempdir()
#     AllEvents_RNASeq_MP<-EventDetectionMultipath(SG_RNASeq,cores=1,Path=TxtPath,paths=3)

## ----ListofLists_MP, eval=FALSE-----------------------------------------------
#  Events[[i]][[j]]

## ----reclassification, eval=TRUE----------------------------------------------
#load splicing graph
data("SG_reclassify")

#load table with info of the events
PathFiles<-system.file("extdata",package="EventPointer")
inputFile <- paste(PathFiles,"/Events_found_class.txt",sep="")
EventTable <- read.delim(file=inputFile)
#this table has the information of 5 complex events.

EventTable_new <- Events_ReClassification(EventTable = EventTable,
                        SplicingGraph = SG_reclassify)

## ----PSI_ADV, eval=TRUE, collapse=TRUE----------------------------------------

# Microarrays (two paths)
data(ArraysData)
PSI_Arrays_list<-EventPointer:::getPSI(ArraysData)
PSI_Arrays <- PSI_Arrays_list$PSI
Residuals_Arrays <- PSI_Arrays_list$Residuals

# Microarrays (Multi-Path)
data(ArrayDatamultipath)
PSI_MP_Arrays_list <- EventPointer:::getPSImultipath(ArrayDatamultipath)
PSI_multipath_Arrays <- PSI_MP_Arrays_list$PSI
Residuals_multipath_Arrays <- PSI_MP_Arrays_list$Residuals

# RNASeq (two paths)
data(AllEvents_RNASeq)
PSI_RNASeq_list<-EventPointer:::getPSI_RNASeq(AllEvents_RNASeq)
PSI_RNASeq <- PSI_RNASeq_list$PSI
Residuals_RNASeq <- PSI_RNASeq_list$Residuals

# RNASeq (Multi-Path)
data(AllEvents_RNASeq_MP)
PSI_MP_RNASeq_list <- EventPointer:::getPSI_RNASeq_MultiPath(AllEvents_RNASeq_MP)
PSI_multipath_RNASeq <- PSI_MP_RNASeq_list$PSI
Residuals_multipath_RNASeq <- PSI_MP_RNASeq_list$Residuals


## ----PSI_ADV2, eval=TRUE, collapse=TRUE---------------------------------------

Dmatrix<-matrix(c(1,1,1,1,0,0,1,1),nrow=4,ncol=2,byrow=FALSE)
Cmatrix<-t(c(0,1))

table <- PSI_Statistic(PSI = PSI_Arrays,Design = Dmatrix,Contrast = Cmatrix,nboot = 20)


## ----PSI_ADV3, eval=TRUE, collapse=TRUE---------------------------------------

Dmatrix<-matrix(c(1,1,1,1,0,0,1,1),nrow=4,ncol=2,byrow=FALSE)
Cmatrix<-t(t(c(0,1)))

Ress <- vector("list", length = ncol(Cmatrix))

fitresiduals <- limma::lmFit(Residuals_Arrays,design = Dmatrix)
fitresiduals2 <- limma::contrasts.fit(fitresiduals, Cmatrix)
fitresiduals2 <- limma::eBayes(fitresiduals2)

# repeated if there is more than one contrast
for (jj in 1:ncol(Cmatrix)) {
  TopPSI <- limma::topTable(fitresiduals2, coef = jj, number = Inf)[, 1, drop = FALSE]
  colnames(TopPSI)<-"Residuals"
  Ress[[jj]] <- TopPSI
}



## -----------------------------------------------------------------------------
sessionInfo()

