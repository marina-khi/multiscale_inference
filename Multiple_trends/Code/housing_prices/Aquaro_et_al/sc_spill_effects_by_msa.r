# (follows from another R script, see below)

#-- compute LeSage's direct, spill-in/out effects --#

# compute partial effects
source("sc_partial_effects.r") #cc#

# loop over explanatory variables
l_effect <- vector("list", length(v_variable)) # (2,)
names(l_effect) <- v_variable
for (variable_shortName in v_variable) {
     #cat(sprintf("  %10s: %s\n", "Variable", variable_shortName))

     a_pe <- l_pe[[variable_shortName]] # (N,N,n_h)

     # loop over h
     l_effect_k <- vector("list", length(v_h)) # (n_h,)
     for (i_h in seq_along(v_h)) {
          h <- v_h[[i_h]]
          # matrix of partial effects
          m_pe <- a_pe[,, i_h]

          # direct effects
          v_direct_effects <- diag(m_pe) # (N,)

          # spill-in/-out effects (cumulative)
          diag(m_pe) <- 0 # handy when computing the sum of off-diagonal elements
          v_cumulative_spill_in_effects <- rowSums(m_pe)
          v_cumulative_spill_ou_effects <- colSums(m_pe)

          # store everything in a tibble
          df_effect_k_h <- tibble::tibble(
               fips       = v_fips, 
               effect_dir = v_direct_effects, 
               effect_csi = v_cumulative_spill_in_effects, 
               effect_cso = v_cumulative_spill_ou_effects,
               h          = h) # (N,5)

          # store
          l_effect_k[[i_h]] <- df_effect_k_h
     }
     # put all dfs into a single df, one on top of the other
     df_effect_k <- bind_rows(l_effect_k) # (NH,5)

     # store
     l_effect[[variable_shortName]] <- df_effect_k
}

# save
saveRDS(l_effect, file = "./data2/effect_spill_by_msa.rds")
