rm(list = ls())

#-- compute partial effects --#

library(dplyr)

cat("\n  ***\n  Calculating partial effects...\n")

v_h <- seq(from = 0, to = 10, by = 1) #cc#
v_variable <- c("pop", "inc") #cc#

# load HSAR estimates
df_theta <- readr::read_csv(file = "./data2/estimates_W075.csv", col_names = TRUE) # (338,?)

# keep columns of interest
df_theta <- df_theta %>% select(fips, region, Wy0, pp, ic, Wy1, y1) # (338,7)
v_fips          <- df_theta[["fips"  ]] # (338,)
v_N_region_code <- df_theta[["region"]] # (338,)

# load W (spatial weight matrix)
df_W <- readr::read_csv("./data1/yang/W75.csv", col_names = FALSE) # (377,377)
m_W <- as.matrix(df_W)

# remove MSAs without neighbours from W
v_rowSums <- rowSums(m_W)
v_ind <- which(v_rowSums == 0)
if (length(v_ind) > 0) {
     m_W <- m_W[-v_ind, ] 
     m_W <- m_W[, -v_ind] 
     cat(sprintf("\n  N = %d\n", nrow(m_W)))
}

# compute partial effects
## retrieve estimates
v_psi0    <- df_theta[["Wy0"]] # (N,)
v_psi1    <- df_theta[["Wy1"]]
v_lambda  <- df_theta[["y1" ]]
v_beta_pp <- df_theta[["pp" ]]
v_beta_in <- df_theta[["ic" ]]

N <- nrow(df_theta)
m_Psi0    <- diag(v_psi0   , nrow = N) # (N,N)
m_Psi1    <- diag(v_psi1   , nrow = N)
m_Lambda  <- diag(v_lambda , nrow = N)
m_Beta_pp <- diag(v_beta_pp, nrow = N)
m_Beta_in <- diag(v_beta_in, nrow = N)

m_Phi <- solve(diag(N) - (m_Psi0 %*% m_W), m_Lambda + (m_Psi1 %*% m_W))

# loop over explanatory variables
l_pe <- vector("list", length(v_variable)) # (2,)
names(l_pe) <- v_variable
for (variable_shortName in v_variable) {
     #cat(sprintf("  %10s: %s\n", "Variable", variable_shortName))

     if (variable_shortName == "pop") {
          m_Beta <- m_Beta_pp
     } else if (variable_shortName == "inc") {
          m_Beta <- m_Beta_in
     }

     # loop over h
     n_h <- length(v_h)
     a_pe <- array(NA_real_, dim = c(N, N, n_h))
     for (i_h in seq_along(v_h)) {
          h <- v_h[[i_h]]

          # matrix power
          m_Phi_h <- matrixcalc::matrix.power(m_Phi, h)

          # matrix of partial effects
          m_pe <- m_Phi_h %*% solve(diag(N) - (m_Psi0 %*% m_W), m_Beta) # (N,N)

          # store
          a_pe[,, i_h] <- m_pe

          # write values to a file for Natalia
          if (0) { #cc#
               if (h == 0) {
                    df_foo <- tibble::tibble(fips = v_fips, region = v_N_region_code) # (338,2)
                    df_foo <- df_foo %>% 
                         mutate(region = if_else(region == 1, 2, region)) %>%
                         mutate(region = if_else(region == 6, 7, region))
                    df_boo <- tibble::as_tibble(m_pe) # (338,338)
                    names(df_boo) <- v_fips
                    df_bar <- bind_cols(df_foo, df_boo) # (338,338+2)
                    readr::write_csv(
                         x = df_bar, 
                         path = sprintf("./data2/partial_effects_%s_h%d.csv", variable_shortName, h))
               }
          }
     }
     # store
     l_pe[[variable_shortName]] <- a_pe
}
