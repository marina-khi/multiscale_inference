Code documentation for 
"Multiscale Inference and Long-Run Variance Estimation in Nonparametric Regression with Time Series Errors"

Authors:
Marina Khismatullina, University of Bonn, marina.k@uni-bonn.de
Michael Vogt, Universiy of Bonn, michael.vogt@uni-bonn.de

=============================================================================================================
Data
=============================================================================================================

Hadley Centre Central England Temperature (HadCET) dataset,
Monthly Mean Central England Temperature (Degrees C) 

The CET dataset is the longest instrumental record of temperature in the world. It contains the mean monthly surface air temperatures (in degrees Celsius) from the year 1659 to the present.
These monthly temperatures are representative of a roughly triangular area of the United Kingdom enclosed by Lancashire, London and Bristol. Manley (1953, 1974) compiled most of the monthly series, covering 1659 to 1973.
These data were updated to 1991 by Parker et al (1992). It is now kept up to date by the Climate Data Monitoring section of the Hadley Centre, Met Office.
Since 1974 the data have been adjusted to allow for urban warming: currently a correction of -0.2 °C is applied to mean temperatures.

CET datasets are freely available for use under Open Government License. They can be downloaded from https://www.metoffice.gov.uk/hadobs/hadcet/.
--------------------------------------------------------------------------------------------------------------

Hadley Centre Climatic Research Unit Temperature (HadCRUT4) dataset,
Global Monthly and Annual Temperature Anomalies (Degrees C) 

HadCRUT4 is a global temperature dataset, providing gridded temperature anomalies across the world. This dataset has been developed by the Climatic Research Unit (University of East Anglia)
in conjunction with the Hadley Centre (UK Met Office). It is being continually updated and expanded by P. Jones of the Climatic Research Unit (CRU).
The temperature series covers the period of 1850-2015 and the anomalies are calculated relative to the 1961-90 reference period means.

HadCRUT4 is subject to Crown copyright protection. The material may be downloaded to file or printer for the purposes of private study and scientific research.
Any other proposed use of the material is subject to a copyright licence available from the Met Office. The data can be downloaded from https://cdiac.ess-dive.lbl.gov/trends/temp/jonescru/jones.html.

=============================================================================================================
Files
=============================================================================================================

The overall structure of the code is as follows. There are four main files each of which produces a specific part of the simulations and applications:

- main_size.r produces the size simulations for our multiscale test and dependent SiZer reported in Section 5.1.1.
- main_power.r produces the power simulations for our multiscale test and dependent SiZer reported in Section 5.1.2.
- main_lrv.r produces the simulation results for our long-run variance estimator, the estimator of Hall and Van Keilegom (2003) and the oracle estimator reported in Section 5.2.
- main_app.r produces the application results from Section 6, where our multiscale test is applied to UK and global temperature data.


These main files read in a number of functions which are collected in the folder .\functions. The simulation and application results are stored either as figures or as .tex files (for tables) in the folder .\plots.
The tables and figures are as in the paper up to seed.

Each main file is divided into several blocks, each block being responsible for one part of the calculations with either a table or a series of plots as a result. Each block has a title that shortly describes
what it is responsible for. The blocks are separated by a series of hashes (#) and are independent of each other. If you want to run only one specific part of the calculations,
you need to run the code in the very beginning of the main file (that contains all the references to the libraries and auxiliary functions) and then the code of the corresponding block. You do not need to run previous blocks.

In order to run the code on your computer, you will need R packages Rcpp and xtable. Theses packages are freely available on CRAN.

All programs are written in R with some functions in C++. They are all quite self-explanatory and commented. The code is self-sufficient and the parts from C++ are read into the R code automatically. 
