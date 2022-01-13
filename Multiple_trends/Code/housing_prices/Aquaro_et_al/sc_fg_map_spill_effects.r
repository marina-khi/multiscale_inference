rm(list=ls())

#-- thematic maps of LeSage's direct, spill-in/out effects --#

library(dplyr)
library(magrittr)
library(ggplot2)
library(sf)

# -------------------------------------------------------------------------- #
#cc#
v_h <- c(0, 3, 6)

# -------------------------------------------------------------------------- #
# load shapefiles
df_sf_msa <- readRDS(file = "./data1/shapes_msa.rds")
df_sf_sta <- readRDS(file = "./data1/shapes_sta.rds")

# -------------------------------------------------------------------------- #
# load direct/spill-in/spill-out effects
l_effect <- readRDS(file = "./data2/effect_spill_by_msa.rds") # (2,)

# -------------------------------------------------------------------------- #
# maps
# loop over explanatory variables
for (variable_shortName in c("pop", "inc")) { # Population, Income
for (lesage_shortName in c("dir", "csi", "cso")) { # Direct, Cumulative spill-in, Cumulative spill-out
     
#for (variable_shortName in c("pop")) { #cc#
#for (lesage_shortName in c("csi")) {

variable_longName <- NULL
lesage_longName   <- NULL
df_effect <- NULL

# labels
if (variable_shortName == "pop") {
     variable_longName <- "Population"
} else if (variable_shortName == "inc") {
     variable_longName <- "Income"
}

# select: pop (inc)
df_effect <- l_effect[[variable_shortName]] # (338*3,??)
df_effect %<>% filter(h %in% v_h)
df_effect[["var_numeric"]] <- NULL

# select: direct (spill-in/-out)
if (lesage_shortName == "dir") {
     lesage_longName <- "Direct"
     df_effect %<>% rename(var_numeric = effect_dir)
} else if (lesage_shortName == "csi") {
     lesage_longName <- "Cumulative spill-in"
     df_effect %<>% rename(var_numeric = effect_csi)
} else if (lesage_shortName == "cso") {
     lesage_longName <- "Cumulative spill-out"
     df_effect %<>% rename(var_numeric = effect_cso)
}
df_effect %<>% select(fips, h, var_numeric)

cat(sprintf("  %10s: %s\n", "Variable", variable_longName))
cat(sprintf("  %10s: %s\n", "Effects" , lesage_longName))

# transform a vector from "numeric" to discrete-"character"
## prepare breaks (bins)
v_quartiles <- quantile(
     abs(df_effect$var_numeric), 
     probs = c(0.25, 0.50, 0.75))

q25 <- v_quartiles[["25%"]]
q50 <- v_quartiles[["50%"]]
q75 <- v_quartiles[["75%"]]

## from continous to discrete with 9 categories (4 below zero, 4 above). The 9th category ("No-neighbours) comes after left_join()
# Note: I'm using "character" instead of "integer" for two reasons:
# 1. Types other than "character" (such as "integer") are coerced to "character" before being transformed into "factor", see later in the code
# 2. I need to add an extra chategory later ("No-Neigh"), and this operation is easier/cleaner with "character" than with "factor" 
v_breaks1 <- c(-Inf, -q75, -q50, -q25, 0, q25, q50, q75, Inf) # (9,)

df_effect %<>% mutate(var_inte = cut(var_numeric, breaks = v_breaks1, labels = FALSE)) %>%
     mutate(var_char = as.character(var_inte))

# check that the categories above were exhaustive
stopifnot(!is.na(df_effect$var_char))

for (h in v_h) {
# select rows corresponding to h == h
df_effects_h <- df_effect[df_effect$h == h, ] # (338,4)

# merge with spatial coordinates
df_sf_effect_h <- left_join(df_sf_msa, df_effects_h, by = c("CBSAFP" = "fips")) # (377,5)
stopifnot(nrow(df_sf_effect_h) == nrow(df_sf_msa))

# transform NAs (here MSA without neighbours) to an extra level in a factor variable
df_sf_effect_h %<>% mutate(var_char = if_else(is.na(var_char), "9", var_char))

# from "character" to "factor"
v_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
df_sf_effect_h %<>% mutate(var_fact = factor(var_char, levels = v_levels, ordered = TRUE))

# -------------------------------------------------------------------------- #
# maps
map <- ggplot()

# US states
map <- map + geom_sf(
     data        = df_sf_sta,
     fill        = NA,
     color       = "dark grey",
     size        = 0.1)

# MSAs
map <- map + geom_sf(
     data    = df_sf_effect_h,
     mapping = aes(fill = var_fact),
     color   = "light grey", 
     size    = 0.1,
     inherit.aes = FALSE)

## scale_discrete_manual()
### named colour vector
v_values <- c(
     "1" = "#B2182B",
     "2" = "#D6604D",
     "3" = "#F4A582",
     "4" = "#FDDBC7",
     "5" = "#D1E5F0",
     "6" = "#92C5DE",
     "7" = "#4393C3",
     "8" = "#2166AC",
     "9" = "gray90")

### labels
v_labels <- vector("character", length(v_values))
v_labels[1] <- sprintf("<=%.3f", -q75)
v_labels[2] <- sprintf("(%.3f,%.3f]", -q75, -q50)
v_labels[3] <- sprintf("(%.3f,%.3f]", -q50, -q25)
v_labels[4] <- sprintf("(%.3f,0]", -q25)
v_labels[5] <- sprintf("(0,%.3f]", q25)
v_labels[6] <- sprintf("(%.3f,%.3f]", q25, q50)
v_labels[7] <- sprintf("(%.3f,%.3f]", q50, q75)
v_labels[8] <- sprintf(">%.3f", q75)
v_labels[9] <- "No-Neigh"

map <- map + scale_discrete_manual(
     aesthetics = c("fill"),
     values     = v_values, # "B2182B"
     breaks     = v_levels, # "1", ?redundant if var is factor?
     labels     = v_labels, # "(,]"
     drop       = FALSE # Should unused factor levels be omitted from the scale?
     )

# house keeping
map <- map + 
     guides(fill = guide_legend(
          nrow           = 1,
          title.position = "top",
          title.hjust    = 0.5,
          label.position = "bottom",
          label.hjust    = 0.5)) +
     theme(
         axis.text         = element_blank(),
         axis.ticks        = element_blank(),
         panel.background  = element_blank(),
         panel.grid        = element_blank(),
         legend.title      = element_text(size = 10), # I've made these a bit larger
         legend.text       = element_text(size = 6),
         legend.key.height = unit(2, "mm"),
         legend.key.width  = unit(12, "mm"),
         legend.direction  = "horizontal",
         legend.position   = "bottom") +
     labs(
        x    = NULL,
        y    = NULL,
        fill = sprintf("%s, %s, h = %d", variable_longName, lesage_longName, h))

# print
#plot(map)

# save
if (1) {
     filename <- sprintf("./tabfig/fg_map_%s_%s_h%d", variable_shortName, lesage_shortName, h)
     ggsave(
          filename = sprintf("%s.png", filename),
          plot     = map,
          device   = "png",
          type     = "cairo",
          dpi      = 300)
     ggsave(
          filename = sprintf("%s.pdf", filename),
          plot     = map,
          device   = cairo_pdf)
# save
}
}}} # end for-loops
