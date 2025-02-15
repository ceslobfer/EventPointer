# EventPointer

[![platforms](http://bioconductor.org/shields/availability/3.8/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)
[![build](http://bioconductor.org/shields/build/release/bioc/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)
[![updated](http://bioconductor.org/shields/lastcommit/release/bioc/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)
[![yearsBioC](http://bioconductor.org/shields/years-in-bioc/EventPointer.svg)](http://bioconductor.org/packages/release/bioc/html/EventPointer.html)

The *EventPointer* R package offers a streamlined method for identifying, classifying, and visualizing alternative splicing events using RNA-seq data. There are two primary workflows: EventPointerBAM (EP_BAM), which utilizes BAM files from splice-aware aligners (e.g., STAR) to detect, classify, and quantify splicing events, allowing for de novo event detection; and EventPointerST (EP_ST), which relies on GTF transcriptome annotations and quantification from pseudoaligners (e.g., salmon, kallisto) for event detection and classification. Each workflow is divided into two steps: detecting, cataloging, and calculating the Percent Spliced In (PSI) of events, and performing statistical analysis to calculate PSI increments between conditions.

![**Figure 1.** EventPointer pipeline ](./vignettes/EPST_EPBAM.png)

# Installation
EventPointer can be installed from Bioconductor using the BiocManager package:

```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

library(BiocManager)

BiocManager::install("EventPointer")
```
