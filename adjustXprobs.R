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
