############################
#Calculating size and power#
############################

results_size     <- simulations_size(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)
results_power    <- simulations_power(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)


###################################
#Implementing clustering algorithm#
###################################

simulations_clustering(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)
results_clusters <- clustering_analysis(N_ts_sim, different_T, different_alpha)