# Figure 1

This folder contains the scripts for reproducing Figure 1 of the manuscript.

- [metadataProcessing.Rmd](metadataProcessing.Rmd) is used to explore the number of epigenomic experiments processed
by the Ensembl Regulation primary analysis pipeline to summarise it for Figure 1a.

- [RFsSummary.Rmd](RFsSummary.Rmd) is used to plot the number of features per species identified by the Ensembl Regulatory Build displayed in Figure 1b.

---

## Requirements

- [R](https://cran.r-project.org/) (>= 4.0) with packages:
  - [dplyr](https://cloud.r-project.org/web/packages/dplyr/index.html)
  - [stringr](https://cloud.r-project.org/web/packages/stringr/index.html)
  - [tidyr](https://cloud.r-project.org/web/packages/tidyr/index.html) 
  - [ggplot2](https://cloud.r-project.org/web/packages/ggplot2/index.html)
  - [cowplot](https://cloud.r-project.org/web/packages/cowplot/index.html)
  - [ggh4x](https://cloud.r-project.org/web/packages/ggh4x/index.html)
