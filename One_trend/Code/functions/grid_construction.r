
grid_construction <- function(T, u.grid=NULL, h.grid=NULL, boundaries=NULL){


  # This function computes the location-bandwidth grid for the multiscale test. 
  #
  # Inputs:
  # T           sample size
  # u.grid      vector of location points in the unit interval [0,1]
  #             If u.grid=NULL, a default grid is used.
  # h.grid      vector of bandwidths
  #             If h.grid=NULL, a default grid is used.
  # boundaries  If "delete", all combinations of locations u and bandwidths h 
  #             for which [u-h,u+h] is not a subset of [0,1] are deleted.
  #             Default is boundaries=NULL. 
  # 
  # Outputs: list with the elements
  # gset        matrix of all location-bandwidth points (u,h) that remain 
  #             after deletion of boundary points (if boundaries = "delete"); 
  #             the i-th row gset[i,] corresponds to the i-th point (u,h) 
  # bws         vector of bandwidths (after deletion of boundary points)
  # lens        vector of length=length(bws), lens[i] gives the number of 
  #             locations in the grid for the i-th bandwidth level
  # grid.type   specifies whether the default grid is used or not
  # gset.full   matrix of all location-bandwidth points (u,h) including 
  #             points at the boundary
  # pos.full    vector indicating which points (u,h) have been deleted  
  #             (only needed to produce plots later on) 
 

  grid.type <- "non-default"
  if(is.null(u.grid) | is.null(h.grid))
  { grid.type <- "default"
    u.grid    <- seq(from = 5/T, to = 1, by = 5/T)
    h.grid    <- seq(from = 3/T, to = 1/4+3/T, by = 5/T)
  }

  gset      <- expand.grid(u=u.grid, h=h.grid) 
  gset.full <- gset
  N         <- dim(gset)[1]
  pos.vec   <- 1:N

  if(!is.null(boundaries))
  { for(i in 1:N)
    {  if(gset[i,1] - gset[i,2] < 0 | gset[i,1] + gset[i,2] > 1) 
         pos.vec[i] <- NA     
    }
  }

  pos.full <- pos.vec
  pos.vec  <- pos.vec[!is.na(pos.vec)]
  gset     <- gset[pos.vec,] 

  bws <- unique(gset[,2])
  lengths.u <- rep(NA,length(bws))
  for(i in 1:length(bws))
     lengths.u[i] <- sum(gset[,2] == bws[i])

  return(list(gset=gset, bws=bws, lens=lengths.u, gtype=grid.type, gset_full=gset.full, pos_full=pos.full))

}
   

