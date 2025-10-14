## Experiments analysis

This folder contains the scripts for reproducing the Figure 1 of the manuscript. 

- [metadataProcessing.Rmd](./scripts/metadataProcessing.Rmd) is used for exploring the number of epigenomic experiments processed
by the Ensembl Regulation [primary analysis pipeline](https://regulation.ensembl.org/114/help/primary_processing) and summarising it for Figure 1a.

- [RFsSummary.Rmd](./scripts/RFsSummary.Rmd) is used to download regulatory features files identified with the Ensembl [Regulatory Build](https://regulation.ensembl.org/help/regulatory_build) from the [Ensembl 114 FTP](https://ftp.ensembl.org/pub/release-114/regulation/) and summarise them for Figure 1b.

---

### Requirements

- [R](https://cran.r-project.org/) (>= 4.0) with packages: 
  - [argparse](https://cloud.r-project.org/web/packages/argparse/index.html) 
  - [dplyr](https://cloud.r-project.org/web/packages/dplyr/index.html)
  - [stringr](https://cloud.r-project.org/web/packages/stringr/index.html)
  - [tidyr](https://cloud.r-project.org/web/packages/tidyr/index.html) 
  - [bedtoolsr](http://phanstiel-lab.med.unc.edu/bedtoolsr.html)
