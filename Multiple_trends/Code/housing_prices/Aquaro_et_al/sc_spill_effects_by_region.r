rm(list = ls())

#-- compute average direct, spill-in/out effects, _by region_ --#

# compute partial effects
source("sc_partial_effects.r") #cc#

# group regions with few MSAs
v_N_region_code[v_N_region_code == 1] <- 2 # NE + ME
v_N_region_code[v_N_region_code == 6] <- 7 # SW + RM

v_R_region_code <- c(2, 3, 4, 5, 7, 8) #cc#
n_R <- length(v_R_region_code)
     
# regional indices
l_ind <- vector("list", n_R)
for (i_r in 1:n_R) {
     r <- v_R_region_code[[i_r]]
     v_ind <- which(v_N_region_code == r) # (N_r,1)
     stopifnot(length(v_ind) > 0)
     l_ind[[i_r]] <- v_ind
}

# compute average partial effect, by region

# loop over pop and inc
l_effect <- vector("list", length(v_variable)) # (2,)
names(l_effect) <- v_variable
for (variable_shortName in v_variable) {
     #cat(sprintf("  %10s: %s\n", "Variable", variable_shortName))

     a_pe <- l_pe[[variable_shortName]] # (N,N,n_h)

     # loop over h
     l_effect_k <- vector("list", length(v_h))
     for (i_h in seq_along(v_h)) {
          h <- v_h[[i_h]]

          # matrix of partial effects
          m_pe <- a_pe[,, i_h]

          # loop over regions
          v_direct1_effect <- vector("numeric", n_R)
          v_direct0_effect <- vector("numeric", n_R)
          m_sel <- matrix(NA_real_, n_R, n_R) # "Sum of ELements" of the sub-matrix
          m_nel <- matrix(NA_real_, n_R, n_R) # "Num of ELements" of the sub-matrix
          for (ri in 1:n_R) {
               v_ind_i <- l_ind[[ri]]
               for (rj in 1:n_R) {
                    v_ind_j <- l_ind[[rj]]

                    # extract sub-matrix corresponding to regions (i,j)
                    m_pe_ij <- m_pe[v_ind_i, v_ind_j] # (N_ri,N_rj)

                    if (ri == rj) {
                         # average direct effect
                         # (average over elements on the main diagonal)
                         v_direct1_effect[[ri]] <- mean(diag(m_pe_ij)) # (1,)
                         
                         # average indirect effect
                         # (average over off diagonal elements)
                         m_pe_ij_2 <- m_pe_ij
                         diag(m_pe_ij_2) <- 0
                         n_row <- nrow(m_pe_ij_2)
                         v_direct0_effect[[ri]] <- sum(m_pe_ij_2) / ((n_row * n_row) - n_row) # (1,)
                    } else {
                         # these quantities are computed only for off-diagonal sub-matrices
                         m_sel[ri, rj] <- sum(m_pe_ij) # (1,)
                         m_nel[ri, rj] <- nrow(m_pe_ij) * ncol(m_pe_ij) # (1,)
                    }
               }
          }

          # average cumulative spill-in/-out effects
          ## handy when computing the sum of off-diagonal elements
          diag(m_sel) <- 0
          diag(m_nel) <- 0

          v_cumulative_spill_in_effect <- rowSums(m_sel) / rowSums(m_nel)
          v_cumulative_spill_ou_effect <- colSums(m_sel) / colSums(m_nel)

          # store everything in a tibble
          df_effect_k_h <- tibble::tibble(
               region_code = v_R_region_code, 
               effect_di1  = v_direct1_effect, 
               effect_di0  = v_direct0_effect, 
               effect_csi  = v_cumulative_spill_in_effect, 
               effect_cso  = v_cumulative_spill_ou_effect,
               h           = h) # (R,6)

          # store
          l_effect_k[[i_h]] <- df_effect_k_h
     }

     # put all dfs into a single df, one on top of the other
     df_effect_k <- bind_rows(l_effect_k) # (HR,6)

     # store
     l_effect[[variable_shortName]] <- df_effect_k
}

# save
saveRDS(l_effect, file = "./data2/effect_spill_by_region.rds")
