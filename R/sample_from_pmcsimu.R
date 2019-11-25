#!/opt/local/bin/Rscript

# Darren Wraith, Martin Kilbinger 2008, 2012
# Resampling from a weighted point distribution (PMC simulation)


library(MASS)
library(getopt)
library(methods)
library(optparse)


parser = OptionParser(usage = "sample_from_pmcsimu.R [options] sample [sample1 [sample2 ...]]") #, option_list = option_list)

parser = add_option(parser, c("-M", "--mult"), type="integer", default=1,
  help="Output sample MULT times input (default %default)")
parser = add_option(parser, c("-i", "--input"), default="pmcsim",
  help="Input PMC simulation file (default '%default')")
parser = add_option(parser, c("-o", "--output"), default="sample",
  help="Output sample file (default '%default')")
parser = add_option(parser, c("-s", "--sample_method"), type="integer", default=0,
  help="(Re)sample method. 0: weighted (default), 1: residual, 2: systematic)")
parser = add_option(parser, c("", "--nclipw"), type="integer", default=0,
  help="Clip points with NCLIWP largest weights (default 0)")


cl = parse_args(parser, positional_arguments = TRUE)


inname   = cl$option$input
outname  = cl$option$output
resample = cl$option$sample_method
mult     = cl$options$mult
nclipw   = cl$options$nclipw


x1 = x2 = NULL;
systema = function(we, N, mult) {
	x1 = trunc(N * mult * cumsum(we) + runif(1))
	x2 = (1:N)[(x1>0) & (!duplicated(x1))]
	rep(x2, c(x1[x2[1]], x1[x2[-1]] - x1[x2[-length(x2)]]))
}


psim    = read.table(inname)
dimf    = dim(psim)[2]
nsample = length(psim[,1])
#prob    = exp(psim[,1])

# MKDEBUG new
pmax   = max(psim[,1])
prob   = exp(psim[,1] - pmax)

# Clip points with highest weights
if (nclipw>0) {
  iclip = order(prob)[(nsample - nclipw + 1) : nsample]
  #cat (iclip, "\n")
  cat ("Setting ", nclipw, " highest weights to zero: ", log(prob[iclip]), "\n")
  prob[iclip] = 0
}



if (resample==0) {
  
  # Usual, weighted resampling
  cat("Weighted resampling (default method)\n")
  indices = sample(1:nsample, size=nsample * mult, prob=prob, replace=T)

} else if (resample==1) {

  # Residual resampling
  cat("Residual resampling\n")
  rw = floor(nsample * prob)
  indices1 = rep.int(1:nsample, times=rw)
  suppressWarnings(sample)
  cat(sum(rw))
  cat("\n")
  cat(nsample)
  cat("\n")
  indices2 = sample(1:nsample, size=(nsample - sum(rw)) * mult, prob=(nsample * prob - rw), replace=T)
  indices  = c(indices1, indices2)
} else if (resample==2) {
  
  # Systematic resampling
  cat("Systematic resampling\n")
  indices = systema(prob, nsample, mult)
  
} else {
  
  cat("Unknown resample type\n")
  stop()
  
}


sample = matrix(NA, nrow=nsample, ncol=dimf)
sample = psim[indices,]
write.table(sample, outname, col.names=F, row.names=F)

# end

