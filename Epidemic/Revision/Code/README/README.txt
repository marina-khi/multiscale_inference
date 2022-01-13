Code documentation for 
"Simultaneous statistical inference for epidemic trends"

Authors:
Marina Khismatullina, University of Bonn, marina.k@uni-bonn.de
Michael Vogt, Universiy of Bonn, michael.vogt@uni-bonn.de

=============================================================================================================
Data
=============================================================================================================

Data on the geographic distribution of COVID-19 cases worldwide (© ECDC [2005-2019])

The dataset is updated daily and contains the latest available public data on COVID-19. It contains the number of new cases (and deaths) reported per day and per country.
The first available date is 31 December 2019.

The dataset can be downloaded from https://www.ecdc.europa.eu/en/publications-data/download-todays-data-geographic-distribution-covid-19-cases-worldwide.

It is owned by ECDC and is subject to ECDC’s copyright policy. The dataset is public and may be reproduced, adapted and/or distributed, totally or in part,
irrespective of the means and/or the formats used, provided that ECDC is always acknowledged as the original source of the material.
Such acknowledgement must be included in each copy of the material. 

--------------------------------------------------------------------------------------------------------------

The Oxford COVID-19 Government Response Tracker (OxCGRT)

OxCGRT is a dataset that contains data on various policy responses to the pandemic of COVID-19. This dataset has been developed by Blavatnik School of Government and is continuously updated
by a cross-disciplinary Oxford University team. As of 9 July, the dataset has data from more than 160 countries. It includes 17 most common indicators such as school closures,
travel restictions and income support for citizens from 1 January 2020 to the present. The data from these indicators is then aggregated into a set of four indices:
an overall government response index, a containment and health index, an economic support index and the original stringency index. A full description of the data and how it is collected,
you can find in the working paper ‘Variation in government response to COVID-19’ that can be accessed freely at https://www.bsg.ox.ac.uk/research/publications/variation-government-responses-covid-19.

OxCGRT dataset is freely available and can be downloaded from https://www.bsg.ox.ac.uk/research/research-projects/coronavirus-government-response-tracker.

=============================================================================================================
Files
=============================================================================================================

The overall structure of the code is as follows. There are four main files each of which produces a specific part of the simulations and applications:

- main_sim.r produces the size and power simulations for our multiscale test reported in Section 4.1.
- main_app.r produces the application results from Section 4.2, where our multiscale test is applied to the data on new daily cases of COVID-19.

Our multiscale test is implemented in the R package multiscale, available from GitHub at ...
In order to run the code on your computer, you will need our R package multiscale as well as R packages Rcpp and xtable. The latter packages are freely available on CRAN.

These main files read in a number of functions which are collected in file .\functions\functions.r. The simulation and application results are stored either as figures or as .tex files (for tables) in the folder .\plots.
The tables and figures are as in the paper up to the seed.

The main file for size and power simulations (main_sim.r) is divided into several blocks, each block being responsible for one part of the calculations with either a table as a result.
Each block has a title that shortly describes what it is responsible for. The blocks are separated by a series of hashes (#) and are independent of each other.
If you want to run only one specific part of the calculations, you need to run the code in the very beginning of the main file (that contains all the references to the libraries and auxiliary functions)
and then the code of the corresponding block. You do not need to run previous blocks.

All programs are written in R. They are all quite self-explanatory and commented. The code is self-sufficient. 
