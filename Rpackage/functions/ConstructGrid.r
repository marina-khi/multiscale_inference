grid_construction <- function(T, u.grid=NULL, h.grid=NULL, deletions=NULL){
  # This function computes the location-bandwidth grid for the multiscale test. 
  #
  # Inputs:
  # T           sample size
  # u.grid      vector of location points in the unit interval [0,1]
  #             If u.grid=NULL, a default grid is used.
  # h.grid      vector of bandwidths
  #             If h.grid=NULL, a default grid is used.
  # deletions   vector of location-bandwidth points (u,h) that are deleted
  #             from the grid. Default is deletions=NULL in which case 
  #             nothing is deleted.
  # 
  # Outputs: list with the elements
  # gset        matrix of location-bandwidth points (u,h) that remains 
  #             after deletions, the i-th row gset[i,] corresponds to the 
  #             i-th point (u,h) 
  # bws         vector of bandwidths (after deletions)
  # lens        vector of length=length(bws), lens[i] gives the number of 
  #             locations in the grid for the i-th bandwidth level
  # gtype       specifies whether the default grid is used or not
  # gset_full   matrix of all location-bandwidth points (u,h) including 
  #             deleted points
  # pos_full    vector indicating which points (u,h) have been deleted  
 

  grid.type <- "non-default"
  if(is.null(u.grid) & is.null(h.grid) & is.null(deletions))
  { grid.type <- "default"
    #Tlen <- min(T,1000)
    Tlen <- T
    u.grid <- seq(from = 1/Tlen, to = 1, by = 1/Tlen)
    h.grid <- seq(from = 1/Tlen, to = 1/4, by = 1/Tlen)
    h.grid <- h.grid[h.grid > log(Tlen)/Tlen]
  }

  gset      <- expand.grid(u=u.grid, h=h.grid) 
  gset.full <- gset
  N         <- dim(gset)[1]
  pos.vec   <- 1:N

  if(!is.null(deletions))
    pos.vec <- deletions

  pos.full <- pos.vec
  pos.vec  <- pos.vec[!is.na(pos.vec)]
  gset     <- gset[pos.vec,] 

  bws <- unique(gset[,2])
  lengths.u <- rep(NA,length(bws))
  for(i in 1:length(bws))
     lengths.u[i] <- sum(gset[,2] == bws[i])
  
  correct = sqrt(2 * log(1 / (2 * gset[, 2])))   
  
  return(list(gset=gset, bws=bws, lens=lengths.u, gtype=grid.type, gset_full=gset.full, pos_full=pos.full, correct = correct))
}