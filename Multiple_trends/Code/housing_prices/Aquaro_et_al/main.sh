#!/usr/bin/env bash

# defatoring
matlab -batch 'sc_defactoring'
# estimating the HSAR model
matlab -batch 'sc_hsar'
# Tab. 3
matlab -batch 'sc_tb_mg_byregion'
# save the results as csv to be used in R
matlab -batch 'sc_mat2csv_estimates_hsar'
# Fig. 1, 2, and F1
Rscript sc_fg_map.r
# generate spill-over effects at MSA- and at regional-level
Rscript sc_spill_effects_by_msa.r
Rscript sc_spill_effects_by_region.r
# Fig. F2-F4
Rscript sc_fg_map_spill_effects.r
# Tab. F1
Rscript sc_tb_density_by_region.r
# Tab. F2
Rscript sc_tb_spill_effects_by_region_dir_ind.r

# -------------------------------------------------------------------------- #
# Tab. F3
# 377x377 distantance matrix
matlab -batch 'sc_distance_bw_MSAs'
# 377x377 inverse distance matrix
matlab -batch 'sc_generate_W_inverse_distance'

# before running these last 3 lines, change "miles = 75" to "miles = 0" in the 3 scripts below
#matlab -batch 'sc_defactoring'
#matlab -batch 'sc_hsar'
#matlab -batch 'sc_tb_mg_byregion'
