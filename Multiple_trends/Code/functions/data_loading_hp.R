##############
#DATA LOADING#
##############

library(readxl)
library(tidyr)
library(dplyr)
library(seasonal)
library(zoo)


hp1 <- read_excel("data/JSTdatasetR5.xlsx", sheet = "Data", na = "") #Macroeconomic dataset for 1870-2017
hp2 <- read_dta("data/NPLH.dta")                                     #Knoll et al. dataset

hp  <- merge(select(hp1, 'year', 'iso', 'country', 'cpi', 'rgdppc', 'pop', 'ltrate'),
             select(hp2, 'year', 'iso', 'hpnom'), by = c('iso', "year"), all = TRUE)


hp_data <-
  hp %>%
  subset((year <= 2012) & (year >= 1890)) %>%
  mutate(hpreal = hpnom / (cpi / 100)) %>%   #Calculating real housing prices
  group_by(iso) %>%
  mutate(infl = (cpi - dplyr::lag(cpi, n = 1, default = NA)) / 100) %>%
  mutate(infl = coalesce(infl, 0)) %>%       #Inflation as changes in CPI
  transform(hpreal = na.approx(hpreal)) %>%  #Imputing the missing values
  transform(ltrate = na.approx(ltrate)) %>%
  mutate(log_hp = log(hpreal)) %>%
  mutate(log_gdp = log(rgdppc)) %>%
  mutate(log_pop = log(pop)) %>%
  group_by(iso) %>%
  mutate(delta_log_hp = log(hpreal) - log(dplyr::lag(hpreal, n = 1, default = NA))) %>%
  mutate(delta_log_pop = log(pop) - log(dplyr::lag(pop, n = 1, default = NA))) %>%
  mutate(delta_log_gdp = log(rgdppc) - log(dplyr::lag(rgdppc, n = 1, default = NA))) %>%
  mutate(delta_ltrate = ltrate - dplyr::lag(ltrate, n = 1, default = NA)) %>%
  mutate(delta_infl = infl - dplyr::lag(infl, n = 1, default = NA)) %>%
  subset(iso %in% c('AUS', 'BEL', 'DNK', 'FRA', 'NLD', 'NOR', 'SWE', 'USA'))

rm(list=setdiff(ls(), c("hp_data", "alpha", "sim_runs", "q", "r")))