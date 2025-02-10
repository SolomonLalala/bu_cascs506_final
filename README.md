# bu_cascs506_final
## Description
Aging is complex biological process with decline in physiology and increased vulnerability to diseases. Developing intervention on aging and aging related diseases involves measuring biological age and aging rate at the molecular level. Aging clocks are machine learning models trained with molecular features such as genes or proteins to estimate an individual's age. Epigenetic aging clock build upon DNA methylation pattern that changes with aging. Epigenetic aging clock is robust in predicting age across various tissue types and can measure epigenetic changes shared between aging and many cancers. 

## Goals
- Predict the aging rate of a cohort based on DNA methylation data. 
- (If time allows, predict cancer risk of a cohort based on DNA methylation data.)

## Data Collection
DNA methylation data will be downloaded from the database below. Will be primarily using Gene Expression Omnibus database.
- **Gene Expression Omnibus**:  [GEO](https://www.ncbi.nlm.nih.gov/geo/)
  - **GSE40279**:  [GSE40279](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279)
  - **GSE55763**:  [GSE55763](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55763)
  - **GSE87571**:  [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571)
  - **GSE167998**:  [GSE167998](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167998)
- **The Cancer Genome Atlas Program**: [TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga)
- **Features**: DNA methylation maker and age.

## Visualization 

## Modelling
Penalized linear regression models

## Test Plan

## Reference
Rutledge J, Oh H, Wyss-Coray T. 2022. Measuring biological age using omics data. Nat Rev Genet. 23(12):715–727. doi:10.1038/s41576-022-00511-7.
Hannum G, Guinney J, Zhao L, Zhang L, Hughes G, Sadda S, Klotzle B, Bibikova M, Fan J-B, Gao Y, et al. 2013. Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. Mol Cell. 49(2):359–367. doi:10.1016/j.molcel.2012.10.016.
