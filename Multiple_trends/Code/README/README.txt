Code documentation for 
"Multiscale Comparison of Nonparametric Trend Curves"

Authors:
Marina Khismatullina, Erasmus University Rotterdam, khismatullina@ese.eur.nl
Michael Vogt, Ulm University, m.vogt@uni-ulm.de

=============================================================================================================
Data
=============================================================================================================

Data used in the analysis of the time series trends in GDP growth in Section 7.1 were collected from multiple sources: Refinitiv Datastream, the OECD.Stat database,
Federal Reserve Economics Data (FRED) and the Barro-Lee Educational Attainment dataset. The data were aggregated and stored in the the file /data/gdp_data.RData.
The description of the variables and the sources they were collected from is provided in the main body of the paper.

--------------------------------------------------------------------------------------------------------------

House price indices

Knoll et al. (2917) presents a novel dataset that covers residential house price indices for 14 advanced economies over the years 1870 to 2012. 
The data are distributed through openICPSR, a public access repository supported by the Inter-university Consortium for Political and Social Research (ICPSR), under an Other License.
The dataset is freely available and was downloaded on 13 January 2022 from https://www.openicpsr.org/openicpsr/project/113055/version/V1/view.

The full citation for the dataset is as follows:

Knoll, Katharina, Schularick, Moritz, and Steger, Thomas. Replication data for: No Price Like Home: Global House Prices, 1870-2012.
Nashville, TN: American Economic Association [publisher], 2017. Ann Arbor, MI: Inter-university Consortium for Political and Social Research [distributor], 2019-10-12.
https://doi.org/10.3886/E113055V1

--------------------------------------------------------------------------------------------------------------

The Jordà-Schularick-Taylor Macrohistory Database

The Jordà-Schularick-Taylor Macrohistory Database is one of the most extensive long-run macro-financial dataset to date. It covers 18 advanced economies since 1870 on an annual basis.
The database comprises 48 real and nominal variables such as bank credit to the non-financial private sector, mortgage lending and long-term returns on housing, equities, bonds and bills.
As stated on the website, the database captures the near-universe of advanced-country macroeconomic and asset price dynamics, covering on average over 90 percent of advanced-economy output
and over 50 percent of world output.

A full description of the data and how it is collected, you can find in the following paper: 
Òscar Jordà, Moritz Schularick, and Alan M. Taylor. 2017. Macrofinancial History and the New Business Cycle Facts. In NBER Macroeconomics Annual 2016, volume 31,     
edited by Martin Eichenbaum and Jonathan A. Parker. Chicago: University  of Chicago Press.

The Jordà-Schularick-Taylor Macrohistory Database is freely available may be used, shared and adapted subject to Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).
For full license see https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode. Here is a human-readable summary: https://creativecommons.org/licenses/by-nc-sa/4.0/.                
 
The data were downloaded on 13 January 2022 from https://www.macrohistory.net/database/.


=============================================================================================================
Files
=============================================================================================================

The overall structure of the code is as follows. There are four main files each of which produces a specific part of the simulations and applications:

- main_sim_test.r produces the size and power simulations for our multiscale test reported in Section 6.
- main_sim_test.r produces the finite sample properties of the clustering algorithm reported in Section 6.
- main_app_gdp.r produces the application results from Section 7.1, where our multiscale test and the clustering procedure are applied to compare the trends in the GDP time series.
- main_app_hp.r produces the application results from Section 7.2, where our multiscale test and the clustering procedure are applied to compare the trends in the real house prices.

Our multiscale test is implemented in the R package "multiscale", available from GitHub at https://github.com/marina-khi/multiscale_inference.
In order to run the code on your computer, you will need our R package "multiscale" as well as R packages "Rcpp", "dplyr", "tidyr", "zoo", "haven", "dendextend", "xtable", "car",
"ggplot2", "Matrix", "foreach", "parallel", "iterators", "doParallel" and "seasonal". The latter packages are freely available on CRAN.

These main files read in a number of functions which are collected in the folder .\functions. The simulation and application results are stored either as figures or as .tex files (for tables)
in the folder .\output and the subfolders therein. The tables and figures are as in the paper up to the seed.

The main files are divided into several blocks, each block being responsible for one part of the calculations. Each block has a title that shortly describes what it is responsible for.
The blocks are separated by a series of hashes (#).

All programs are written in R. They are all quite self-explanatory and commented. The code is self-sufficient. 
