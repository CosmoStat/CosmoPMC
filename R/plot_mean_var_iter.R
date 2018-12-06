# MK 6/2008
# Plots mean and errors (1,2,3sigma) for a PMC or MCMC simulation.
# MCMC: give nsamples, fsfinal, niter on the command line

library(Hmisc, warn.conflicts = F)


if (file.exists("config_pmc")==T) {
   configname = "config_pmc"
} else {
   if (file.exists("config_mcmc")==T) {
      configname = "config_mcmc"
    }
}


# Get parameters from config file
niter = nsamples = fsfinal = -1
config <- readLines(configname)
for (ln in 1:length(config)) {

	if (pmatch("npar", config[ln], nomatch=0)!=0) {
		tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
		npar = as.double(tmp[2])
	}
	if (pmatch("niter", config[ln], nomatch=0)!=0) {
		tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
		niter = as.double(tmp[2])
	}
	if (pmatch("min", config[ln], nomatch=0)!=0) {
		tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
		min = as.double(tmp[2:(npar+1)])
	}

	if (pmatch("max", config[ln], nomatch=0)!=0) {
		tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
		max = as.double(tmp[2:(npar+1)])
	}
	if (pmatch("nsamples", config[ln], nomatch=0)!=0) {
		tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
		nsamples = as.double(tmp[2])
	}
	if (pmatch("fsfinal", config[ln], nomatch=0)!=0) {
		tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
		fsfinal = as.double(tmp[2])
	}

}

# If no config_pmc, use command  line arguments
args = commandArgs(TRUE)
if (nsamples==-1) { nsamples = as.double(args[1]) }
if (fsfinal==-1) { fsfinal = as.double(args[2]) }
if (niter==-1) { niter = as.integer(args[3]) }

output = system(paste("get_spar.pl -c ", configname, " R"), intern=T)
lab    = unlist(strsplit(output, "&"))

# Read files "mean_<iter>" into 3d-array m
nmom = 7
m = array(1:niter*npar*nmom, c(niter,npar,nmom))
for (iter in 1:niter) {
	filename = paste("mean", iter-1, sep="_")
	tmp = read.table(filename)
	for (par in 1:npar) {
		for (mom in 1:nmom) {
			m[iter,par,mom] = tmp[par,mom+1]
		}
	}
}


# Plot
#x = (1:niter)*nsamples
#x[niter] = x[niter]*fsfinal
x = (1:niter)
x[niter] = x[niter]*fsfinal

for (par in 1:npar) {

	tmp = paste("mean_var_iter", par-1, sep="_");
	outname = paste(tmp, "ps", sep=".")
	postscript(outname, encoding = "TeXtext.enc", paper="special", width=5, height=5, horizontal=F)

	mean = m[,par,1]
	plot(x, mean, ylim=c(min[par],max[par]), xlab=paste("sample points (", expression(10^4), ")", sep=""),
	     ylab=parse(text=lab[par]))

	# 1sigma
	sig1p = m[,par,2]
	sig1m = m[,par,3]
	for (iter in 1:niter) {
		if (sig1p[iter]<0) { sig1p[iter] = max[par]-mean[iter] }
		if (sig1m[iter]<0) { sig1m[iter] = mean[iter]-min[par] }
	}

	# 2sigma
	sig2p = m[,par,4]
	sig2m = m[,par,5]
	for (iter in 1:niter) {
		if (sig2p[iter]<0) { sig2p[iter] = max[par]-mean[iter] }
		if (sig2m[iter]<0) { sig2m[iter] = mean[iter]-min[par] }
	}

	# 3sigma
	sig3p = m[,par,6]
	sig3m = m[,par,7]
	for (iter in 1:niter) {
		if (sig3p[iter]<0) { sig3p[iter] = max[par]-mean[iter] }
		if (sig3m[iter]<0) { sig3m[iter] = mean[iter]-min[par] }
	}

	errbar(x, mean, mean+sig1p, mean-sig1m, add=T)
	#errbar(x, mean, mean+sig2p, mean-sig2m, add=T, lty=2)
	errbar(x, mean, mean+sig3p, mean-sig3m, add=T, lty=3)

	dev.off()

}

# end
