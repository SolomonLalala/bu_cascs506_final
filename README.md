# bu_cascs506_final
###### Aging is complex biological process with decline in physiology and increased vulnerability to diseases. Developing intervention on aging and aging related diseases involves measuring biological age and aging rate at the molecular level. Aging clocks are machine learning models trained with molecular features such as genes or proteins to estimate an individual's age. Epigenetic aging clock build upon DNA methylation pattern that changes with aging. Epigenetic aging clock is robust in predicting age across various tissue types and can measure epigenetic changes shared between aging and many cancers. This project uses existing methods developed by Hannum et al. to train an epigenetic aging clock and aims to predict the age and aging rate of individuals using public available DNA methylation data.

# Final report
- For presentation, see <https://youtu.be/HIhTsgAyCFs> I set this to be private and shareed with all staff. Invitees need sign in to their Google Account (bu email) to view the private video.
- For report, view final_report.html or [Click here to view the final report html](https://SolomonLalala.github.io/bu_cascs506_final/final_report.html)
- To reproduce the result
1.  Install R (4.4.3) and Rstudio. 
2.  Clone the repository. 
- For the full pipeline, open final_report.Rmd and run the codes.
- For testing, I tried to build a GitHub workflow to run the test.R but failed. Therefore, for testing code, please run: Rscript -e 'install.packages(c("tidyverse", "knitr", "caret", "glmnet", "ggplot2", "GEOquery", "foreach", "doParallel", "testthat")) in bash and then run test.R.

# Reference

###### Anastasiadi D, Piferrer F. 2023. Bioinformatic analysis for age prediction using epigenetic clocks: Application to fisheries management and conservation biology. Front Mar Sci. 10. <doi:10.3389/fmars.2023.1096909>.

###### Farrell CP. A Simple Epigenetic Clock Using Python and SciKit-Learn. Colin P Farrell. <https://colinpfarrell.com/tag-ec-tutorial/>.

###### Hannum G, Guinney J, Zhao L, Zhang L, Hughes G, Sadda S, Klotzle B, Bibikova M, Fan J-B, Gao Y, et al. 2013. Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates. Mol Cell. 49(2):359–367. <doi:10.1016/j.molcel.2012.10.016>.

###### OpenAI. (2025). ChatGPT (Mar 2025 version). Accessed via <https://chat.openai.com>
