
# compute default version of location-bandwidth grid 

default_grid <- function(T)
{   u.pts <- seq(from = 5/T, to = 1, by = 5/T)
    h.pts <- seq(from = 3/T, to = 1/4+3/T, by = 5/T)
    # full grid
    # u.pts <- seq(from = 1/T, to = 1, by = 1/T)
    # h.pts <- seq(from = 3/T, to = 1/4+1/T, by = 1/T)
    grid <- expand.grid(u = u.pts, h = h.pts)  
    # subsetting u and h such that [u-h, u+h] lies in [0,1]
    # grid <- subset(grid, u >= 0 & u <= 1, select = c(u, h)) 
    return(grid) 
}