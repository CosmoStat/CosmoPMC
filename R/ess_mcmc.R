# Martin Kilbinger, Darren Wraith 2008
# Effective sample size (ess) for an MCMC chain

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

n = length(chain[1,]) - 2
sample = chain[,3:n]

ess = effectiveSize(mcmc(sample))
cat(paste(ess, " "))
cat("\n")
q()
