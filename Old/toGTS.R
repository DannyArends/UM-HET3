

toGT <- function(gtsprob, marker){
  mp <- gtsprob[, grep(marker, colnames(gtsprob))]
  gts <- unlist(lapply(lapply(lapply(apply(mp, 1, function(x){which(x > 0.85)}),names), strsplit, ":"), function(x){
    if(length(x) > 0){ return(x[[1]][2]); }else{ return(NA) }
  }))
  return(gts)
}

a <- toGT(gtsp, "1_3010274")
b <- toGT(gtsp, "c1.loc1")
x <- cbind(a[which(!a == b)],  b[which(!a == b)])

a <- toGT(gtsp, "X_36008085")
b <- toGT(gtsp, "cX.loc1")
x <- cbind(a[which(!a == b)],  b[which(!a == b)])

plot(gtsp[, "X_36008085:BC"] - gtsp[, "cX.loc1:BC"])
plot(gtsp[, "1_3010274:BC"] - gtsp[, "c1.loc1:BC"])
