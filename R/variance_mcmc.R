# Darren Wraith, Martin Kilbinger 2008
# Variance of the mean of an MCM chain

# Argument: chain


library(lattice)
library(coda)

args = commandArgs(TRUE)

name = args[1]
if (file.exists(name)) {
  chain = read.table(name)
} else {
  cat(paste("File ", name, " not found"))
  q()
}

d  = dim(chain)[2]
n1 = length(as.ts(chain[,1]))

for (i in 3:d) {
  p0      = spectrum0(as.ts(chain[,i]))
  varmean = p0$spec[1]/n1
  str     = sprintf("%3d %.4e\n", i-3, varmean)
  cat(str)
}

q()
