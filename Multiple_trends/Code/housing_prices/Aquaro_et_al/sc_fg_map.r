rm(list = ls())

#-- Maps with MSA-specific estimates --#

library(dplyr)
library(magrittr)
library(ggplot2)
library(sf)

# -------------------------------------------------------------------------- #
#cc#
verbose <- 1

# -------------------------------------------------------------------------- #
# load shapefiles
df_sf_msa <- readRDS(file = "./data1/shapes_msa.rds")
df_sf_sta <- readRDS(file = "./data1/shapes_sta.rds")

# -------------------------------------------------------------------------- #
# load HSAR estimates
df_theta <- readr::read_csv(file = "./data2/estimates_W075.csv", col_names = TRUE)

# keep variables of interest
df_theta %<>% select(fips, region, Wy0, pp, ic, Wy1, y1, Wy0_Wy1)

# check that Wy0+Wy1 is in (-1,1) otherwise breaks need to be adjusted accordingly #cc#
if (max(abs(df_theta[["Wy0_Wy1"]])) >= 1) {
     error("  Adjust breaks accordingly")
}

# boundary cases
df_theta %<>% mutate(boundary_case = (abs(Wy0) > 0.994) + (abs(y1) > 0.994) + (abs(Wy1) > 0.994)) # 0, 1, 2, 3
if (verbose) {
     # number of boundary cases
     ans <- nrow(df_theta %>% filter(boundary_case > 0))
     cat(sprintf("\n  Number of |psi_i| > 0.994 OR |lambda_i| > 0.994 OR |psi_1i| > 0.994: %d\n", ans))
     
     # table with regional MG estimates
     cat(sprintf("\n  Regional mean-group estimates (excluding boudary cases):\n"))
     df_ans <- df_theta %>% 
          filter(boundary_case == 0) %>% 
          select(-c(fips, boundary_case)) %>%
          mutate(region = if_else(region == 1, 2, region)) %>%
          mutate(region = if_else(region == 6, 7, region)) %>%
          group_by(region) %>% 
          select(region, Wy0_Wy1, Wy0, Wy1, y1, pp, ic) %>%
          summarise_all(mean) %>%
          mutate_if(is.double, list(~round(., 3)))
     print(df_ans)
     
     # table with US MG estimates
     cat(sprintf("\n  US mean-group estimates (excluding boudary cases):\n"))
     df_ans <- df_theta %>% 
          filter(boundary_case == 0) %>% 
          select(-c(fips, region, boundary_case)) %>% 
          select(Wy0_Wy1, Wy0, Wy1, y1, pp, ic) %>%
          summarise_all(mean) %>%
          mutate_if(is.double, list(~round(., 3)))
     print(df_ans)
}

# -------------------------------------------------------------------------- #
# loop
for (variable in c("Wy0", "pp", "ic", "Wy1", "y1", "Wy0_Wy1")) {
     s_label_name_variable  <- NULL
     v_breaks               <- NULL
     v_labels               <- NULL
     df_theta[["var_nume"]] <- NULL
     df_theta[["var_inte"]] <- NULL
     df_theta[["var_char"]] <- NULL
     df_theta[["var_fact"]] <- NULL

# -------------------------------------------------------------------------- #
# select variable of interest
if (variable == "Wy0") {
     df_theta[["var_nume"]] <- df_theta[["Wy0"]]
     s_label_name_variable <- "$\\hat{\\psi}_{0i}$"
} else if (variable == "pp") {
     df_theta[["var_nume"]] <- df_theta[["pp"]]
     s_label_name_variable <- "$\\hat{\\beta}_{i}^{\\pop}$"
} else if (variable == "ic") {
     df_theta[["var_nume"]] <- df_theta[["ic"]]
     s_label_name_variable <- "$\\hat{\\beta}_{i}^{\\inc}$"
} else if (variable == "y1") {
     df_theta[["var_nume"]] <- df_theta[["y1"]]
     s_label_name_variable <- "$\\hat{\\lambda}_{i}$"
} else if (variable == "Wy1") {
     df_theta[["var_nume"]] <- df_theta[["Wy1"]]
     s_label_name_variable <- "$\\hat{\\psi}_{1i}$"
} else if (variable == "Wy0_Wy1") {
     df_theta[["var_nume"]] <- df_theta[["Wy0_Wy1"]]
     s_label_name_variable <- "$\\hat{\\psi}_{0i}+\\hat{\\psi}_{1i}$"
}

# -------------------------------------------------------------------------- #
# breaks
# prepare breaks (bins) and labels to discretise the estimates
# from continous to discrete with 10 categories:
# 4 below zero and 4 above (with variable-specific breaks)
# 1 "boundary cases"
# 1 "No-neighbours", which comes after left_join()

if (variable == "pp") {
     v_breaks <- c(-Inf, seq(from = -1.5, to = 1.5, by = 0.5), Inf)
     v_labels <- c("<-1.5", "(-1.5,-1]", "(-1,-0.5]", "(-0.5,0]", "(0,0.5]", "(0.5,1]", "(1,1.5]", ">1.5")
} else if (variable == "ic") {
     v_breaks <- c(-Inf, c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), Inf)
     v_labels <- c("<-0.3", "(-0.3,-0.2]", "(-0.2,-0.1]", "(-0.1,0]", "(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", ">0.3")
} else {
     v_breaks <- seq(from = -1, to = 1, by = 0.25)
     v_labels <- c("(-1,-0.75]", "(-0.75,-0.5]", "(-0.5,-0.25]", "(-0.25,0]", "(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1)")
}
stopifnot(length(v_breaks) == 9) # 8+1
stopifnot(length(v_labels) == 8)
v_labels <- c(v_labels, "Non-Conv", "No-Neigh")

df_theta %<>% mutate(var_inte = cut(var_nume, breaks = v_breaks, labels = FALSE)) %>% # from continuous to discrete
     mutate(var_char = as.character(var_inte)) %>% # from "integer" (e.g. 1) to "character" (e.g. "1")
     mutate(var_char = if_else(boundary_case > 0, "9", var_char)) # boundary cases

# -------------------------------------------------------------------------- #
# merge with spatial coordinates
df_sf_theta <- left_join(df_sf_msa, df_theta, by = c("CBSAFP" = "fips"))
stopifnot(nrow(df_sf_theta) == nrow(df_sf_msa))

# transform NAs (here MSA without neighbours) to an extra level in a factor variable
df_sf_theta %<>% mutate(var_char = if_else(is.na(var_char), "10", var_char)) # No-Neigh

# final touch: from "character" to "factor"
v_levels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10") # this is more clear than as.character(seq(.))
#df_sf_theta %<>% mutate(var_fact = factor(var_char, levels = v_levels, labels = v_labels, ordered = TRUE)) # 2020-4-13: not sure why label= create probles with scale_discrete_manual()
df_sf_theta %<>% mutate(var_fact = factor(var_char, levels = v_levels, ordered = TRUE))

if (verbose) {
     cat(sprintf("\n  Variable: %s\n", variable))
     cat("  Number of estimates per category:\n")
     print(table(df_sf_theta[["var_fact"]]))

     # glimpse
     #df <- df_sf_theta %>% st_set_geometry(NULL) %>% select(CBSAFP, region, boundary_case, starts_with("var_"))
     #tableHTML::tableHTML(df)
}

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
     data    = df_sf_theta,
     mapping = aes(fill = var_fact),
     color   = "light grey", 
     size    = 0.1,
     inherit.aes = FALSE)

# scale_fill_manual()
## named colour vector
v_values <- c(
     "1" = "#B2182B",
     "2" = "#D6604D",
     "3" = "#F4A582",
     "4" = "#FDDBC7",
     "5" = "#D1E5F0",
     "6" = "#92C5DE",
     "7" = "#4393C3",
     "8" = "#2166AC",
     "9" = "gray80",
     "10" = "gray90")

## as with other scales, "breaks" can be used to control the appearance of the legend.
map <- map + scale_discrete_manual(
     aesthetics = "fill",
     values = v_values, 
     breaks = v_levels,
     labels = v_labels,
     drop   = FALSE) # Should unused factor levels be omitted from the scale?

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
         legend.title      = element_text(size = 8),
         legend.text       = element_text(size = 4),
         legend.key.height = unit(2, "mm"),
         legend.key.width  = unit(10, "mm"),
         legend.direction  = "horizontal",
         legend.position   = "bottom"
         ) +
     labs(
        x    = NULL,
        y    = NULL,
        fill = latex2exp::TeX(s_label_name_variable)
        )

# print
#plot(map)

# save
if (1) {
     filename <- sprintf('./tabfig/fg_map_%s', variable)
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
}
} # end loop
