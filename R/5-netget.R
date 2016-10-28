netget <- function(mat, d, nid, tid, time, 
                   get.graph = FALSE, drop.iso = TRUE, thresh = NULL) {
  obs <- d[, tid] %in% time
  new.d <- d[obs, ]
  s <- min(which(obs))
  e <- max(which(obs))
  m <- mat[s:e, s:e]
  if (!is.null(thresh)) {
    cp <- quantile(m[which(m > 0)], thresh)
    m[which(m < cp)] <- 0
  }
  dimnames(m) <- list(as.character(d[obs, nid]), 
                      as.character(d[obs, nid]))
  if (get.graph) {
    g <- graph.adjacency(m, mode = 'directed', weighted = TRUE, diag = FALSE)    
    if(drop.iso)  {
      drop.v <- which(degree(g) < 1)
      new.d <- new.d[-drop.v, ]
      g <- delete.vertices(g, drop.v)
    }
  }
  else g <- NULL
  list(graph = g, matrix = m, data = new.d)
}