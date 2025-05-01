# bu_cascs506_final

## Final report

[Click here to view the final report](https://SolomonLalala.github.io/bu_cascs506_final/final_report.html)

## Midterm report

midterm report video record: <https://youtu.be/AajjHc2R3Zk> I set this to be private and shareed with all staff. Invitees need sign in to their Google Account (bu email) to view the private video.

## Description

###### Aging is complex biological process with decline in physiology and increased vulnerability to diseases. Developing intervention on aging and aging related diseases involves measuring biological age and aging rate at the molecular level. Aging clocks are machine learning models trained with molecular features such as genes or proteins to estimate an individual's age. Epigenetic aging clock build upon DNA methylation pattern that changes with aging. Epigenetic aging clock is robust in predicting age across various tissue types and can measure epigenetic changes shared between aging and many cancers. This project uses existing methods developed by Hannum et al. to train an epigenetic aging clock and aims to predict the age and aging rate of individuals using public available DNA methylation data.

## Goals
- ###### Predict the age and aging rate of individuals based on DNA methylation data.

## Data Collection

###### DNA methylation data can be downloaded from from Gene Expression Omnibus [GEO](https://www.ncbi.nlm.nih.gov/geo/).

- ### **Primary dataset for training and testing the model**:
    - #### [**GSE40279**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279): Download Series Matrix File using `GEOquery`

- ### **Additional datasets or databases considered for testing model**:

    - #### [**GSE87571**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571)

    - #### [**GSE55763**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55763)

    - #### [**GSE167998**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167998)

    - #### **The Cancer Genome Atlas Program ([TCGA](https://www.cancer.gov/ccg/research/genome-sequencing/tcga))**

- ### **Features**:

  - #### DNA methylation (beta values): continuous variables between 0 and 1, representing the methylation level of a CpG site

    - #### Age: 19 to 101

    - #### Covariates (sex, ethnicity, body mass index, etc)

## Data Visualization

- ### A histogram of the age distribution for samples

- ### A heat map of the age-related DNA methylation across age

## Data Processing

- ### Parse relevant metadata (age, sex, ethnicity, etc) and DNA methylation matrix from the series matrix file

- ### Merge metadata with DNA methylation matrix

- ### Stratify training/testing split (75/25) using age and gender

## Modelling

- ### Apply Elastic Net regression, a penalized regression that have those DNA methylation markers contributing less to prediction excluded by having their coefficients to or close to zero, thus reduce the number of variable

- ### Use 10-fold crossvalidation to select optimal lambda, the regularization parameter

- ### Perform bootstrap 500 times to select DNA methylation markersrobustly associated with age in the dataset

- ### Train final model using optimal lambda and selected DNA methylation markers

## Testing

- ### 656 samples from [**GSE40279**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279) are splited into training (75%) and testing (25%) set by stratified random sampling to preserve the distribution of age and sex

- ### Apply final model on the testing set

- ### Calculate R-squared and RMSE to evaluate the model performance

## Priliminary Results

- ### An elastic net regression model was trained with lambda = 0.4867 and 23 DNA methylation markers. The model predicted the age of individuals with an R-squared of 0.827 and RMSE of 6.04 years.

## Reference

###### Anastasiadi D, Piferrer F. 2023. Bioinformatic analysis for age prediction using epigenetic clocks: Application to fisheries management and conservation biology. Front Mar Sci. 10. <doi:10.3389/fmars.2023.1096909>.

###### Farrell CP. A Simple Epigenetic Clock Using Python and SciKit-Learn. Colin P Farrell. <https://colinpfarrell.com/tag-ec-tutorial/>.

###### Hannum G, Guinney J, Zhao L, Zhang L, Hughes G, Sadda S, Klotzle B, Bibikova M, Fan J-B, Gao Y, et al. 2013. Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. Mol Cell. 49(2):359–367. <doi:10.1016/j.molcel.2012.10.016>.

###### OpenAI. (2025). ChatGPT (Mar 2025 version). Accessed via <https://chat.openai.com>

###### Rutledge J, Oh H, Wyss-Coray T. 2022. Measuring biological age using omics data. Nat Rev Genet. 23(12):715–727. <doi:10.1038/s41576-022-00511-7>.
