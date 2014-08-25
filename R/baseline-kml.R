runKml <- function(data,t.range,n.clusters,n.bins=10,n.restarts=10) {
  breaks <- seq(t.range[1],t.range[2],length.out=n.bins+1)
  expandSequence <- function(s) {
    t <- as.integer(cut(times(s),breaks))
    y <- values(s)
    y.out <- numeric(n.bins)
    for (i in 1:n.bins) {
      y.out[i] <- mean(y[t==i])
    }
    y.out
  }

  seq.ids <- vapply(data,'[[',integer(1),'id')
  exp.seqs <- lapply(data,expandSequence)
  seq.mat <- do.call('rbind',exp.seqs)
  cld <- clusterLongData(seq.mat,seq.ids,1:n.bins,maxNA=floor(n.bins/2))
  invisible(capture.output(kml(cld,n.clusters,n.restarts)))
  cld
}

kmlPartition <- function(cld,k,rank=1) {
  getClusters(cld,k,rank,asInteger=TRUE)
}
