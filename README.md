# Dissection of floral transition by single meristem transcriptomes at high temporal resolution

## Summary
This repository contains code used for key steps in analysis performed in Meir, Aviezer et al. (Nature Plants, 2021), including data processing and generation of figures. 

## Compilation

The SMT_analysis.r script should run independently. However:

  - at each new session, start with sourcing the main functions script:

```r
source("SMT_nat_plants_2021_functions.r")
```

  - Also, make sure you run the scripts from the same parent directory, maintaining folders of processed table (sub-directories Data/ and Supplementary_Tables/) just below it.

  - Some of the analysis rely on files/tables you need to download. The functions that require downloading should put the files in the proper subdirectories. 


A constantly-updating webpage where you can explore the SMT database and look at gene expression kinetics during tomato floral transition:
https://tanaylab.weizmann.ac.il/SMT/

Vignette exemplifying how to order meristems (or other objects) by using “slanter”:
https://tanaylab.github.io/slanter/articles/meristems.html

"slanter" source code and R package:
https://cran.r-project.org/web/packages/slanter/index.html

All raw and processed transcriptome data of this manuscript are deposited in GEO (accession GSE166929):
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166929


