#
# adjustXprobs.R
#
# copyright (c) 2020-2030 - Danny Arends
# 
# Adjust X-chromosome probablitities to avoid issues caused by wrong sex assignment 
#

adjustXprobs <- function(cross){
  sex <- getsex(cross)$sex
  pr <- cross$geno[["X"]]$prob
  stopifnot(!is.null(pr), !is.null(sex))

  for(i in 1:ncol(pr)) {
      pr[sex==0,i,3:4] <- 0
      pr[sex==1,i,1:2] <- 0
      pr[,i,] <- pr[,i,]/rowSums(pr[,i,])
  }
  cross$geno[["X"]]$prob <- pr
  invisible(cross)
}

mlog <- function(msg, append = TRUE){
  cat(paste0(Sys.time(), " - ", msg), file = "log.out", append = append)
  cat(paste0(Sys.time(), " - ", msg))
}
