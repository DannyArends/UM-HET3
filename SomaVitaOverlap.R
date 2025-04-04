## Size of all Mass: 1067 Mb = 44.7
## Size of all Soma: 1298 Mb = 38.1
## Size of all Vita: 872 Mb = 33.5

cnts <- c()
for (x in 1:10000) {
  soma <- sample(1:2600, 30)
  vita <- sample(1:2600, 26)
  cnt <- 0
  for (s in soma) {
    cnt <- cnt + length(which(abs(vita - s) < 10))
  }
  cnts <- c(cnts, cnt)
}

# We find overlapping, 8 Soma peaks within 10 mb of Vita peaks [12 times 95% CI overlap]
(10000 - min(which(sort(cnts) == 8))) / 10000

cnts <- c()
for(x in 1:10000){
  soma <- sample(1:2600, 30)
  mass <- sample(1:2600, 28)
  cnt <- 0
  for(s in mass){
    cnt <- cnt + length(which(abs(vita - s) < 10))
  }
  cnts <- c(cnts, cnt)
}

# We find overlapping, 7 Soma peaks within 10 mb of mass peaks
(10000 - min(which(sort(cnts) == 7))) / 10000
