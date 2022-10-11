rm(list = ls())

#-- calculate population density by region --#

library(dplyr)
library(magrittr)
library(sf)

# -------------------------------------------------------------------------- #
# load Cynthia's data
df_NT <- readr::read_csv("./data1/yang/data_main.csv"
  , col_names = TRUE
  , col_types = readr::cols(
  quarter    = readr::col_character(),
  msacode    = readr::col_integer(),
  regioncode = readr::col_integer(),
  cpi        = readr::col_skip(),
  hp         = readr::col_skip(),
  pop        = readr::col_double(),
  pipc       = readr::col_skip())
  )

# recode regions according to BEA classification
df_NT %<>%
  mutate(regioncode_temp = NA_integer_) %>%
  mutate(regioncode_temp = if_else(regioncode == 4, 1L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 3, 2L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 2, 3L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 5, 4L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 7, 5L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 8, 6L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 6, 7L, regioncode_temp)) %>%
  mutate(regioncode_temp = if_else(regioncode == 1, 8L, regioncode_temp)) %>%
  select(-regioncode) %>%
  rename(regioncode = regioncode_temp)

df_N <- df_NT %>%
     group_by(msacode) %>%
     summarise(
          pop        = median(pop),
          regioncode = first(regioncode)) # (377,3)

# -------------------------------------------------------------------------- #
# load shapefiles (CSBSA)
df_sf_cbsa <- read_sf(dsn = "./data1/shapefiles/tl_2013_us_cbsa.shp", quiet = TRUE)

df_sf_msa <- df_sf_cbsa %>% 
     filter(MEMI == "1") %>% # keep only MSAs
     select(CBSAFP, NAME, ALAND) %>%
     mutate(CBSAFP = as.integer(CBSAFP))

# remove Puerto Rico (PR), Alaska (AK) and Hawaii (HI)
v_id_pr <- c(10380, 11640, 25020, 32420, 38660, 41900, 41980)
v_id_ak <- c(11260, 21820)
v_id_hi <- c(27980, 46520)

v_id_pr_ak_hi <- c(v_id_pr,
                   v_id_ak,
                   v_id_hi)

df_sf_msa %<>% filter(!(CBSAFP %in% v_id_pr_ak_hi)) # (377,4)

# remove geometry
df_msa <- st_set_geometry(df_sf_msa, NULL) # (377,3)

# ALAND: Land area - an area measurement providing the size, in square meters, of the
# land portions of geographic entities for which the Census Bureau tabulates
# and disseminates data.
#
# Area is calculated from the specific boundary recorded for each entity in the
# Census Bureau’s geographic database (see "MAF/TIGER Database").

# Land area measurements are originally recorded as whole square meters (to
# convert square meters to square kilometers, divide by 1,000,000; to convert
# square kilometers to square miles, divide by 2.58999; to convert square
# meters to square miles, divide by 2,589,988).

df_msa %<>% rename(land_square_metres = ALAND) %>%
     mutate(land_square_km = land_square_metres / 1000000) %>%
     mutate(land_square_mi = land_square_km / 2.58999) # (377,6)

# -------------------------------------------------------------------------- #
# put land area and population in the same data-frame
df_merged_N <- left_join(df_N, df_msa, by = c("msacode" = "CBSAFP"))
stopifnot(nrow(df_merged_N) == 377) #cc#

# -------------------------------------------------------------------------- #
# load HSAR estimates
df_theta <- readr::read_csv(file = "./data2/estimates_W075.csv", col_names = TRUE) # (338,?)

# keep variables of interest
df_theta %<>% select(fips) # (338,1)

df_merged_N <- semi_join(df_merged_N, df_theta, by = c("msacode" = "fips")) # (338,1)
     
# -------------------------------------------------------------------------- #
df_merged_N %<>% 
     mutate(regioncode = if_else(regioncode == 1, 2L, regioncode)) %>%
     mutate(regioncode = if_else(regioncode == 6, 7L, regioncode))

# -------------------------------------------------------------------------- #
cat(sprintf("  Number of MSAs: %d\n", nrow(df_merged_N)))
# from MSAs to regions, summing up
df_merged_R <- df_merged_N %>% 
     group_by(regioncode) %>%
     summarize(
          pop = sum(pop),
          n = n(),
          land_square_mi = sum(land_square_mi)) %>%
     mutate(persons_x_square_mi = pop / land_square_mi) # (6,5)

#!! MAKE SURE DATA IS SORTED BEFORE THIS COMMAND
v_R_region_name <- c(
     "New England \\& Mideast",
     "Great Lakes",
     "Plains",
     "Southeast",
     "Southwest \\& Rocky Mountain",
     "Far West")

v_R_region_code2 <- c("1 \\& 2", "3", "4", "5", "6 \\& 7", "8")

df_merged_R %<>% arrange(regioncode) %>%
     mutate(
          regionname  = v_R_region_name,
          regioncode2 = v_R_region_code2) %>%
     select(
          regioncode,
          regioncode2,
          regionname,
          n,
          pop,
          land_square_mi,
          persons_x_square_mi)

readr::write_csv(df_merged_R, path = sprintf("./data2/density_by_region_N%d.csv", nrow(df_merged_N))) 

# write results to a file
closeAllConnections()
fullfilename <- sprintf("./tabfig/in_tb_density_by_region_N%d.tex", nrow(df_merged_N))
con <- file(description = fullfilename, open = "wt")
sink(con, split = TRUE)
for (i_row in 1:nrow(df_merged_R)) {
     for (i_col in 2:ncol(df_merged_R)) {
          ans <- df_merged_R[[i_col]][[i_row]]
          if (i_col == 2) {
               cat(sprintf("%12s", ans))
          } else if (i_col == 3) {
               cat(sprintf(" & %30s", ans))
          } else if (i_col == 4) {
               cat(sprintf(" & %4d", ans))
          } else if (i_col == 5) {
               cat(sprintf(" & %9.0f", ans))
          } else if (i_col == 6) {
               cat(sprintf(" & %7.0f", ans))
          } else if (i_col == 7) {
               cat(sprintf(" & %6.1f", ans))
          }
     }
     cat(" \\\\\n")
}
sink()
close(con)
