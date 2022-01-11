rm(list = ls())

#-- table with average direct, indirect, csi, and cso effects --#

library(dplyr)

v_R_region_name <- c(
"New England \\& Mideast",
"Great Lakes",
"Plains",
"Southeast",
"Southwest \\& Rocky Mountain",
"Far West")

v_R_region_code2 <- c("1 \\& 2", "3", "4", "5", "6 \\& 7", "8")

l_effect <- readRDS(file = "./data2/effect_spill_by_region.rds") # (2,)

v_variable <- names(l_effect)
for (variable_shortName in v_variable) {
     df_effect_k <- l_effect[[variable_shortName]] # (RH,?)

     df_wide <- df_effect_k %>% 
          filter(h %in% c(0, 3, 6)) %>% 
          mutate(effect_di0_percentage = effect_di0 / effect_di1 * 100) %>%
          # select variable of interest
          select(region_code, effect_di1, effect_di0_percentage, h) %>%
          # long to wide (h = 0,3,6)
          tidyr::pivot_wider(
               names_from = h, 
               values_from = c(effect_di1, effect_di0_percentage)
          ) %>%
          arrange(region_code) #%>% select(-region_code)

     # write to a file
     closeAllConnections()
     con <- file(description = sprintf("./tabfig/in_tb_spill_effects_by_region_dir_ind_%s.tex", variable_shortName), open = "wt")
     sink(con, split = TRUE)
     for (i in 1:nrow(df_wide)) {
          cat(sprintf("%6s & %30s", v_R_region_code2[[i]], v_R_region_name[[i]]))
          for (j in 2:ncol(df_wide)) {
               cat(sprintf(" & % 6.3f", df_wide[[j]][[i]]))
          }
          cat(" \\\\\n")
     }
     sink()
}
close(con)
