##############
#DATA LOADING#
##############

library(tidyr)
library(dplyr)
library(seasonal)
library(zoo)

countries <- c("AUT", "AUS", "CAN", "CHE", "DEU", "FIN", "FRA", "GBR",
               "JPN", "NOR", "USA") 

################################
#Loading the human capital data#
################################

h      <- read.csv("data/human_capital.csv", sep = ",", dec = ".",
                   stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")
h$date <- as.Date(paste0("01-01-", h$year), format = "%d-%m-%Y")
h      <- subset(h, select = - c(BLcode, sex, agefrom, ageto, year, region_code))

h_data <- complete(h, date = seq.Date(min(date), max(date), by='quarter'),
                   country)
h_data <- 
  h_data %>%
  arrange(country, date) %>%
  group_by(country) %>%
  fill(WBcode, .direction = "up") %>%
  transform(h_it = na.approx(yr_sch)) %>%
  select(date, WBcode, h_it) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y"))


#########################
#Loading the labour data#
#########################

l           <- read.csv("data/Employment_FRED.csv", sep = ";", dec = ",",
                        stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")
colnames(l) <- c("date", "AUT", "AUS", "BEL", "CAN", "CHE", "DEU",
                 "DNK", "ESP", "FIN", "FRA", "GBR", "GRC",
                 "IRL", "ISL", "ITA", "JPN", "LUX",
                 "NLD", "NOR", "NZL", "PRT", "SWE", "TUR", "USA")
l_data <- 
  l %>%
  mutate(date = as.Date(date, format = "%Y-%m-%d")) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y")) %>%
  complete(date = seq.Date(min(date), max(date), by = 'quarter')) %>%
  arrange(date)

l_data <- gather(l_data, WBcode, l_it_FRED, AUT:USA, factor_key=FALSE)
l_data$l_it_FRED <- l_data$l_it_FRED / 1000 #In thousands of people

l2 <- read.csv("data/datastream_employment.csv", sep = ";", dec = ",",
               stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")

#There are four observations missing for Spain so we fill them in
l2 <-
  l2 %>%
  fill(ESP, .direction = 'up') %>%
  mutate(date = as.Date(as.yearqtr(date, format = 'Q%q %Y'), format = "%d-%m-%Y")) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y"))

l2_data <- gather(l2, WBcode, l_it_datastream, USA:BEL, factor_key=FALSE)

#FINLAND, NORWAY and ITALY - not seasonally adjusted data!!
l3      <- read.csv("data/datastream_employment_not_sa.csv", sep = ";", dec = ",",
                    stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")
l3$date <- as.Date(as.yearqtr(l3$date, format = 'Q%q %Y'), format = "%d-%m-%Y")

#Seasonal adjustment using seasonal package, default options
fin_ts    <- ts(l3$FIN, frequency = 4, start = c(1959, 1))
fin_ts_sa <- seas(x = fin_ts)
l3[!is.na(l3$FIN), 'FIN'] <- fin_ts_sa$data[, 'seasonaladj']

nor_ts    <- ts(l3$NOR, frequency = 4, start = c(1959, 1))
nor_ts_sa <- seas(x = nor_ts)
l3[!is.na(l3$NOR), 'NOR'] <- nor_ts_sa$data[, 'seasonaladj']

ita_ts    <- ts(l3$ITA, frequency = 4, start = c(1959, 1))
ita_ts_sa <- seas(x = ita_ts)
l3[!is.na(l3$ITA), 'ITA'] <- ita_ts_sa$data[, 'seasonaladj']

#There are five observations missing for Italy so we fill them in
l3 <-
  l3 %>%
  fill(ITA, .direction = 'up') %>%
  mutate(date = as.Date(as.yearqtr(date, format = 'Q%q %Y'), format = "%d-%m-%Y")) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y"))

l3_data <- gather(l3, WBcode, l_it_datastream_not_sa, FIN:ITA, factor_key=FALSE)

#PANAMA - monthly data!!
l4 <- read.csv("data/datastream_employment_monthly_sa.csv", sep = ";", dec = ",",
               stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")

l4_data <-
  l4 %>%
  mutate(date = as.Date(date, format = "%d-%m-%Y")) %>%
  group_by(date = paste(quarters(date), lubridate::year(date))) %>%
  summarise(PAN = mean(PAN)) %>%
  mutate(date = as.Date(as.yearqtr(date, format = 'Q%q %Y'), format = "%d-%m-%Y")) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y")) %>%
  arrange(date)

l4_data <- gather(l4_data, WBcode, l_it_datastream_monthly_sa, PAN, factor_key=FALSE)


X_mat <- merge(h_data, l_data, all = TRUE)
X_mat <- merge(X_mat, l2_data, all = TRUE)
X_mat <- merge(X_mat, l3_data, all = TRUE)
X_mat <- merge(X_mat, l4_data, all = TRUE)

X_mat$l_it <- ifelse(X_mat$WBcode %in% c("CHE", "FRA", "BEL", "ESP", "CYP", "CHN", "ISR"),
                     X_mat$l_it_datastream,
                     ifelse(X_mat$WBcode %in% c("FIN", "NOR", "ITA"),
                            X_mat$l_it_datastream_not_sa,
                            ifelse(X_mat$WBcode == 'PAN',
                                   X_mat$l_it_datastream_monthly_sa,
                                   X_mat$l_it_FRED)))

######################
#Loading the gdp data#
######################


gdp <- read.csv("data/gdp_oecd.csv", sep = ",", dec = ".",
                stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")

gdp_data <-
  gdp %>%
  subset(Subject == 'Gross domestic product - expenditure approach') %>%
  mutate(date = as.Date(as.yearqtr(TIME, format = '%Y-Q%q'), format = "%d-%m-%Y")) %>%
  select(date, LOCATION, Value) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y")) %>%
  subset(LOCATION %in% c("AUT", "AUS", "CAN", "CHE", "CHN", "CYP", "DEU",
                         "ESP", "FIN", "FRA", "GBR", "IND", "ISR", "ITA", "JPN",
                         "NOR", "THA", "USA")) %>%
  complete(date = seq.Date(min(date), max(date), by = 'quarter'), LOCATION) %>%
  arrange(LOCATION, date)

colnames(gdp_data) <- c('date', 'WBcode', 'gdp_it_OECD')


gdp_annual           <- read.delim("data/GDP_non_standard_Annual.txt",
                                   stringsAsFactors = FALSE, na.strings = "")
colnames(gdp_annual) <- c("date", "CHN_gdp", "CYP_gdp", "PAN_gdp", "CHN_pop",
                          "CYP_pop", "PAN_pop") 

gdp_annual_data <-
  gdp_annual %>%
  mutate(CHN = CHN_gdp * CHN_pop / 1000000) %>%
  mutate(CYP = CYP_gdp * CYP_pop / 1000000) %>%
  mutate(PAN = PAN_gdp * PAN_pop / 1000000) %>%
  mutate(date = as.Date(date, format = "%Y-%d-%m")) %>%
  select(date, CHN, CYP, PAN) %>%
  subset(date > as.Date("30-09-1974", format = "%d-%m-%Y")) %>%
  complete(date = seq.Date(min(date), max(date), by = 'quarter')) %>%
  arrange(date) %>%
  mutate(CHN = na.approx(CHN) / 4) %>%
  mutate(CYP = na.approx(CYP) / 4) %>%
  mutate(PAN = na.approx(PAN) / 4) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y")) %>%
  arrange(date)

gdp_annual_data <- gather(gdp_annual_data, WBcode, gdp_it_annual, CHN:PAN, factor_key=FALSE)

X_mat <- merge(X_mat, gdp_data, all = TRUE)
X_mat <- merge(X_mat, gdp_annual_data, all = TRUE)

X_mat$gdp_it <- ifelse(X_mat$WBcode %in% c("CHN", "CYP", "PAN"),
                       X_mat$gdp_it_annual, X_mat$gdp_it_OECD)


################################
#Loading the capital stock data#
################################

k           <- read.delim("data/Capital_stock_Annual.txt",
                          stringsAsFactors = FALSE, na.strings = "")
colnames(k) <- c("date", "AUT", "AUS", "BEL", "CAN", "CHE", "DEU",
                 "DNK", "ESP", "FIN", "FRA", "GBR", "GRC",
                 "IRL", "ISL", "ITA", "JPN", "LUX",
                 "NLD", "NOR", "NZL", "PRT", "SWE", "TUR", "USA")
k$date      <- as.Date(k$date, format = "%Y-%m-%d")


k_china           <- read.csv("data/capital_stock_china.csv", sep = ",", dec = ".",
                              stringsAsFactors = FALSE, na.strings = "")
colnames(k_china) <- c("date", "CHN")
k_china$date      <- as.Date(k_china$date, format = "%Y-%m-%d")

k_cyprus           <- read.csv("data/capital_stock_cyprus.csv", sep = ",", dec = ".",
                               stringsAsFactors = FALSE, na.strings = "")
colnames(k_cyprus) <- c("date", "CYP")
k_cyprus$date      <- as.Date(k_cyprus$date, format = "%Y-%m-%d")

k_panama           <- read.csv("data/capital_stock_panama.csv", sep = ",", dec = ".",
                               stringsAsFactors = FALSE, na.strings = "")
colnames(k_panama) <- c("date", "PAN")
k_panama$date      <- as.Date(k_panama$date, format = "%Y-%m-%d")


k <- merge(k, k_china, all = TRUE)
k <- merge(k, k_cyprus, all = TRUE)
k <- merge(k, k_panama, all = TRUE)

k_data <- gather(k, WBcode, Value, AUT:PAN, factor_key = FALSE)

k_data <-
  k_data %>%
  complete(date = seq.Date(min(date), max(date), by='quarter'), WBcode) %>%
  arrange(WBcode, date) %>%
  group_by(WBcode) %>%
  transform(stock_it = na.approx(Value)) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y")) %>%
  select(date, WBcode, stock_it)

X_mat  <- merge(X_mat, k_data, all = TRUE)

#######################
#Loading the gfcf data#
#######################

kf <- read.csv("data/gfcf_oecd.csv", sep = ",", dec = ".",
               stringsAsFactors = FALSE, na.strings = "", fileEncoding="UTF-8-BOM")

kf_data <-
  kf %>%
  mutate(date = as.Date(as.yearqtr(TIME, format = '%Y-Q%q'), format = "%d-%m-%Y")) %>%
  select(date, LOCATION, Value) %>%
  subset(date > as.Date("30-09-1975", format = "%d-%m-%Y")) %>%
  subset(LOCATION %in% c("AUT", "AUS", "CAN", "CHE", "CHN", "CYP", "DEU",
                         "ESP", "FIN", "FRA", "GBR", "IND", "ISR", "ITA", "JPN",
                         "NOR", "THA", "USA")) %>%
  complete(date = seq.Date(min(date), max(date), by = 'quarter'), LOCATION) %>%
  arrange(LOCATION, date)

colnames(kf_data) <- c('date', 'WBcode', 'gfcf_it')

X_mat <- merge(X_mat, kf_data, all = TRUE)

#Some data manipulations
X_mat_filled <-  
  X_mat %>%
  subset(date <= as.Date("30-09-2010", format = "%d-%m-%Y")) %>%
  subset(WBcode %in% countries) %>%
  group_by(WBcode)  %>%
  fill(h_it, .direction = 'down')  #Extrapolating educational attainment for the last two quarters 


#Either we are using capital stock ("stock") or
#gross fixed capital formation ("gfcf") as the variable k:

#X_mat_filled$k_it <- X_mat_filled$stock_it
X_mat_filled$k_it <- X_mat_filled$gfcf_it

save(X_mat_filled, countries, file = "data/gdp_data.RData")
