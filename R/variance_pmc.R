# Darren Wraith, Martin Kilbinger 2008
# Variance of the mean of a PMC simulation

# Argument: pmc_simu

args = commandArgs(TRUE)

name = args[1]
if (file.exists(name)) {
  simul = read.table(name)
} else {
  cat(paste("File ", simul, " not found"))
  q()
}

d  = dim(simul)[2]

for (i in 3:d) {
  varmean = sum(simul[,1]^2*(simul[,i]-sum(simul[,1]*simul[,i]))^2)
  str     = sprintf("%3d %.4e\n", i-3, varmean)
  cat(str)
}

q()
