# Load COSA package
platform="linux"
cosadir="/home/jessez/R/COSA"
source("/home/jessez/R/COSA/r_cosa.q")

# Load data
args <- commandArgs(trailingOnly = TRUE)
filename = args[1]

x <- as.matrix(read.table(filename,header=FALSE))

d <- cosadist(x)
h <- hierclust(d,denplot=F)
labels <- cutree(h,2)

# print labels
for (entry in labels) { 
  if (entry == 1) {print(0)}
  else {print(1)}
}
