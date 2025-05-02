# bu_cascs506_final
###### Aging is complex biological process with decline in physiology and increased vulnerability to diseases. Developing intervention on aging and aging related diseases involves measuring biological age and aging rate at the molecular level. Aging clocks are machine learning models trained with molecular features such as genes or proteins to estimate an individual's age. Epigenetic aging clock build upon DNA methylation pattern that changes with aging. Epigenetic aging clock is robust in predicting age across various tissue types and can measure epigenetic changes shared between aging and many cancers. This project uses existing methods developed by Hannum et al. to train an epigenetic aging clock and aims to predict the age and aging rate of individuals using public available DNA methylation data.

# Final report
- For report, view final_report.html or [Click here to view the final report html](https://SolomonLalala.github.io/bu_cascs506_final/final_report.html)
- To reproduce the result
1.  Install R (4.4.3) and Rstudio. 
2.  Clone the repository. 
- For the full pipeline, open final_report.Rmd and run the codes.
- For testing, I tried to build a GitHub workflow to run the test.R but failed. Therefore, for testing code, please run: Rscript -e 'install.packages(c("tidyverse", "knitr", "caret", "glmnet", "ggplot2", "GEOquery", "foreach", "doParallel", "testthat")) in bash and then run test.R.
