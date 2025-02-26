cnts <- c()
for(x in 1:10000){
  soma <- sample(1:2600, 27)
  vita <- sample(1:2600, 26)
  cnt <- 0
  for(s in soma){
    cnt <- cnt + length(which(abs(vita -s) < 10))
  }
  cnts <- c(cnts, cnt)
}

TOD: REDO Supplment Table 2
