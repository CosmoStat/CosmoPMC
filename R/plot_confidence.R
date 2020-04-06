# Martin Kilbinger, Darren Wraith 2008
# 2013
# Plots 1D density and 2D Confidence levels from a list of sample/chain files

# Arguments: sample1 [sample2 [sample3 [...]]]
# The config file is in read from current directory. If the sample files have a different
# number of parameters, the config file should correspond to the one with the smallest npar.


library(lattice)
library(MASS)
library(coda)
library(getopt)
#library(methods)
library(optparse)
library(funr)



# Command line options
parser = OptionParser(usage = "plot_confidence.R [options]") #, option_list = option_list)

parser = add_option(parser, c("-W","--weights"), action="store_true", default=FALSE,
  help="Use weights in 0th column (default: FALSE, input file is already resampled")
parser = add_option(parser, c("-N", "--Ngrid"), type="integer", default=100,
  help="Number of grid points for smoothing (kde2d) (default %default). Use <=30 for\n\t\tfast-but-dirty plots")
parser = add_option(parser, c("-g", "--gsmooth"), default="30",
  help="Smoothing kernel width, with respect to box size (default %default). In case of more\n\t\tthan one sample, use list separated with '_' for more than value")
parser = add_option(parser, c("-S", "--solid"), action="store_true", default=FALSE,
  help="All contours with solid lines")
parser = add_option(parser, c("-w", "--width"), default="1",
  help="Line width (default 1)")
parser = add_option(parser, c("-k", "--with_keys"), action="store_true", default=FALSE,
  help="Add key to plots")
parser = add_option(parser, c("-K","--keystring"), default="",
  help="Key strings (separate items with '_')")
parser = add_option(parser, c("-L", "--no_key_line"), action="store_true", default=FALSE,
  help="Do not add a line to the keys in the legend")
parser = add_option(parser, c("-c", "--config"), default="config_pmc",
  help="Config file (default '%default')")
parser = add_option(parser, c("-t", "--title"), default="",
  help="Title string for each panel (default empty)")
parser = add_option(parser, c("-b", "--brace"),  action="store_true", default=FALSE,
  help="Add <> as <keystring> (default empty)")
parser = add_option(parser, c("-i", "--index_i"), type="integer", default=-1,
  help="Only create plots with i-th parameter on x-axis")
parser = add_option(parser, c("-j", "--index_j"), type="integer", default=-1,
  help="Only create plots with j-th parameter on y-axis")
parser = add_option(parser, c("-s", "--sigma"), type="integer", default=3,
  help="Plot SIGMA confidence levels (default 3)")
parser = add_option(parser, c("-F", "--color_scheme"), type="integer", default=0,
  help="Color scheme (0 ... 7; default 0)")
parser= add_option(parser, c("-o", "--output"), default="pdf",
 help="Output file format, FORMAT=ps|pdf (default: pdf)")
parser = add_option(parser, c("-m", "--marker"), default="",
  help="Plots a mark at position PAR (e.g. best-fit). Separate items with '_'.")
parser = add_option(parser, c("-e", "--marker_error"), default="",
  help="Plots error bars with rms MARKER_ERROR around the mark (option '-m'). Separate items with '_'.")
parser = add_option(parser, c("--pmin"), default="",
  help="Lower limits for plot. Default: Read from config file (= sample limits). Separate items with '_'.")
parser = add_option(parser, c("--pmax"), default="",
  help="Upper limits for plot. Default: Read from config file (= sample limits). Separate items with '_'.")
parser = add_option(parser, c("-C", "--shade"), default="1",
  help="With shade 0 or 1, default=1")

cl     = parse_args(parser, positional_arguments = TRUE)


use_weights = cl$options$weights
N          = cl$options$N
ssmooth    = cl$options$gsmooth
fsmooth    = as.numeric(unlist(strsplit(ssmooth, "_")))
solid      = cl$options$solid
lwd        = cl$options$width
with_keys  = cl$options$with_keys
keystring  = cl$options$keystring
no_key_line = cl$options$no_key_line
configname = cl$options$config
title      = cl$options$title
index_i    = cl$options$index_i
index_j    = cl$options$index_j
sigma      = cl$options$sigma
color_scheme = cl$options$color_scheme
output_format = cl$options$output
marker     = cl$options$marker
marker_error = cl$options$marker_error
brace      = cl$options$brace
pmins      = cl$options$pmin
pmaxs      = cl$options$pmax
shade     = cl$options$shade

path <- get_script_path()

#tmpname   = "tmptmp.ps"
tmpname   = paste("tmptmp", output_format, sep=".")

blue      =  8/12
green     =  4/12
red       =  0/12
cyan      =  6/12
magenta   = 10/12
yellow    =  2/12
orange    =  3/12

if (color_scheme == 0) {
  mycolors  = c(blue, green, red, rgb(0, 191, 255, maxColorValue=255), rgb(255, 140, 0, maxColorValue=255), magenta)
  strcolors = c("blue", "darkgreen", "red", "deepskyblue", "orange2", "magenta")
} else if (color_scheme == 1) {
  mycolors  = c(green, red, rgb(0, 191, 255, maxColorValue=255), cyan, magenta, blue)
  strcolors = c("darkgreen", "red", "deepskyblue", "cyan", "magenta", "blue")
} else if (color_scheme == 2) {
  mycolors  = c(blue, green, magenta, rgb(255, 140, 0, maxColorValue=255), rgb(0, 191, 255, maxColorValue=255), magenta)
  strcolors = c("blue", "darkgreen", "magenta", "orange2", "deepskyblue", "magenta")
} else if (color_scheme == 3) {
  mycolors  = c(blue, green, magenta, rgb(0, 0, 0, maxColorValue=255), rgb(0, 191, 255, maxColorValue=255), magenta)
  strcolors = c("blue", "darkgreen", "magenta", "black", "deepskyblue", "magenta")
} else if (color_scheme == 4) {
  mycolors  = c(green, magenta, rgb(255, 140, 0, maxColorValue=255), rgb(0, 0, 0, maxColorValue=255), magenta, blue)
  strcolors = c("darkgreen", "magenta", "orange2", "black", "magenta", "blue")
} else if (color_scheme == 5) {
  mycolors  = c(blue, magenta, rgb(0, 0, 0, maxColorValue=255), rgb(0, 191, 255, maxColorValue=255), rgb(0, 191, 255, maxColorValue=255))
  strcolors = c("blue", "magenta", "black", "deepskyblue", "orange2")
} else if (color_scheme == 6) {
  mycolors  = c(orange, blue)
  strcolors = c("orange", "blue")
} else if (color_scheme == 7) {
  mycolors  = c(magenta, blue, green, rgb(0, 0, 0, maxColorValue=255), rgb(0, 191, 255, maxColorValue=255), cyan)
  strcolors = c("magenta", "blue", "darkgreen", "black", "orange2", "cyan")
} else if (color_scheme == 8) {
  mycolors  = c(blue, green, rgb(0, 0, 0, maxColorValue=255), rgb(0, 191, 255, maxColorValue=255), cyan, magenta)
  strcolors = c("blue", "darkgreen", "black", "orange2", "cyan", "magenta")
} else if (color_scheme == 9) {
  mycolors  = c(magenta,  green, rgb(0, 191, 255, maxColorValue=255), rgb(0, 0, 0, maxColorValue=255), blue, cyan)
  strcolors = c("magenta", "darkgreen", "orange2", "black", "blue",  "cyan")
}
prob = c(0.6827)
if (sigma > 1) prob = c(prob, 0.9545)
if (sigma > 2) prob = c(prob, 0.9973)
if (sigma > 3) stop("Sigma (option '-s') too large, has to be 1, 2 or 3")
if (sigma < 1) stop("Sigma (option '-s') too small, has to be 1, 2 or 3")

# Order of colors in HPDregionsplot:
# black, red, green, blue, cyan, magenta, yellow, grey


MK.colors = function (n, alpha = 1, ncol = 1) 
{
  # alpha: transparency
  
  if ((n <- as.integer(n[1])) > 0) {
    k <- n%/%2
    #h <- c(4/12, 2/12, 0/12)     # hue, color
    #s <- c(1, 1, 0)              # saturation, intensity
    #v <- c(0.65, 0.9, 0.95)      # value, lightness
    h <- array(mycolors[ncol], dim=3)
    s <- c(0.0, 0.25, 0.5)
    v <- c(1.0, 0.9, 0.8)
    c(hsv(h = seq.int(h[1], h[2], length.out = k),
          s = seq.int(s[1], s[2], length.out = k),
          v = seq.int(v[1], v[2], length.out = k), 
          alpha = alpha),
      hsv(h = seq.int(h[2], h[3], length.out = n - k + 1)[-1],
          s = seq.int(s[2], s[3], length.out = n - k + 1)[-1],
          v = seq.int(v[2], v[3], length.out = n - k + 1)[-1],
          alpha = alpha))
  }
  else character(0)
}

get_weights = function (psample)
{
  # See sample_from_pmcsimu.R
  pmax   = max(psample[,1])
  prob   = exp(psample[,1] - pmax)
  prob / sum(prob)
}

kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}

# Performs a kernel-density smoothing and plots the density as image (shade=T) or contour (shade=F)
HPDregionplotDW = function (x, vars = 1:2, h = c(1, 1), n = 50, lump = TRUE, prob = 0.95, weights = NULL,
  xlab = NULL, ylab = NULL, ncol = 1, Xlim = NULL, Ylim = NULL, pXlim = NULL, pYlim = NULL, first = T, shade = F, cex.lab=2,
  cex.axis=1.5, ...)
{
  require("MASS")
  parnames <- if (class(x) == "mcmc.list")
    colnames(x[[1]])
  else colnames(x)
  if (is.character(vars)) {
    vars <- match(vars, parnames)
    if (any(is.na(vars)))
      stop("some variable names do not match")
  }
  varnames <- parnames[vars]
  mult <- (class(x) == "mcmc.list" && !lump)
  if (mult)
    stop("multiple chains without lumping not yet implemented")
  if (class(x) == "mcmc.list") {
    if (lump)
      var1 <- c(sapply(x, function(z) z[, vars[1]]))
    else var1 <- lapply(x, function(z) z[, vars[1]])
  }
  else var1 <- x[, vars[1]]
  if (class(x) == "mcmc.list") {
    if (lump)
      var2 <- c(sapply(x, function(z) z[, vars[2]]))
    else var2 <- lapply(x, function(z) z[, vars[2]])
  } else {
    var2 <- x[, vars[2]]
  }


  if (!mult) {
    if (! is.null(weights)) {
      post1 = kde2d.weighted(var1, var2, n = n, h = h, lims = c(Xlim[1], Xlim[2], Ylim[1], Ylim[2]), weights)
    } else {
      post1 = kde2d(var1, var2, n = n, h = h, lims = c(Xlim[1], Xlim[2], Ylim[1], Ylim[2]))
    }
  } else {
    #if (is.list(weights) & length(weights) == 0) {
    if (! is.null(weights)) {
      post1 = mapply(kde2d.weighted, var1, var2, MoreArgs = list(n = n), weights)
    } else {
      post1 = mapply(kde2d, var1, var2, MoreArgs = list(n = n))
    }
  }

  dx = diff(post1$x[1:2])
  dy = diff(post1$y[1:2])
  sz = sort(post1$z)
  c1 = cumsum(sz) * dx * dy
  levels = sapply(prob, function(x) {
    approx(c1, sz, xout = 1 - x)$y
  })
  if (is.null(xlab))
    xlab <- varnames[1]
  if (is.null(ylab))
    ylab <- varnames[2] 

  # For first = TRUE new plot has to be created, Add = F
  if (first==T) {
    Add = F
  } else {
    Add = T
  }

  # Increase margins (if labels get cut off)
  #Default: par(mar=c(5,4,4,2) + 0.1)
  par(mar=c(5,4.25,4,2) + 0.1)

  # Shading out to 99.7%-region
  if (shade==1) {

    prob3 = c(0.6827, 0.9545, 0.9973)
    levels = sapply(prob3, function(x) {
      approx(c1, sz, xout = 1 - x)$y
    })
    Min = levels[length(prob)] # was [3]
    Max = range(post1$z)[2]

    image(post1$x, post1$y, post1$z, c(Min,Max), pXlim, pYlim, col=MK.colors(100, alpha = 1, ncol = ncol),
          xlab = xlab, ylab = ylab, add = Add, cex.lab=cex.lab, cex.axis=cex.axis, las=1)

  } else {
    contour(post1$x, post1$y, post1$z, level = levels, xlab = xlab,
            ylab = ylab, drawlabels = FALSE, add = Add, xlim = pXlim, ylim = pYlim, cex.lab=cex.lab,
            cex.axis=cex.axis, col=strcolors[ncol], las=1, ...)

  }

  invisible(contourLines(post1$x, post1$y, post1$z, level = levels))

  if (Add == FALSE) {
    
    if (pXlim[2] - pXlim[1] < 2) {
      axis(side = 1, at = seq(-10, 10, by = 0.05), labels = FALSE, tcl = -0.2)
    } else {
      axis(side = 1, at = seq(-10, 10, by = 0.2), labels = FALSE, tcl = -0.2)
    }
    
    if (pYlim[2] - pYlim[1] < 2) {
      axis(side = 2, at = seq(-10, 10, by = 0.05), labels = FALSE, tcl = -0.2)
    } else {
      axis(side = 2, at = seq(-10, 10, by = 0.2), labels = FALSE, tcl = -0.2)
    }

  }

}

# Adds the legend to the plot
do_legend = function (keystring, xMin, xMax, yMin, yMax, nsamples, no_line) {
  key_a = unlist(strsplit(keystring, "_"))
  dx    = xMax - xMin
  dxk   = 0.05 * dx
  dy    = yMax - yMin

  dxl   = 0.15 * dx
  dyl   = 0.075 * dy

  fxl   = 0.75 * dx
  fyl   = 0.9  * dy # Top right
  #fyl   = 0.26  * dy # Bottom right

  if (no_line == FALSE) {
    # With lines
    xl = xMin + fxl
    x <- c(xl, xl + dxl)
    off = 0.5
  } else {
    # No lines
    xl = xMin + fxl + dxl + dxk
    off = -1.0
  }
  yl    = yMin + fyl

  for (ns in 1:nsamples) {
    if (no_line == FALSE) {
       y <- c(yl, yl) - (ns-1) * dyl
       if (solid == TRUE) { lty = 1 }
       else { lty = ns }
       lines(x, y, type='l', lty=lty, col=strcolors[ns])
    }

#    key_tmp = expression(""<Map^"2      ", "">"")
#    tmp = '""<Map^"2      ", "">""'
#    key_tmp = eval(expression(tmp))


    if (brace == TRUE) {
      keytmp<- list(bquote(paste(""<"", .(key_a[ns]), "">"", sep="")))
      text(xl - dxk, yl - (ns-1) * dyl,  do.call(expression, keytmp), offset=off, pos='2', col=strcolors[ns])
    } else {
        text(xl - dxk, yl - (ns-1) * dyl, key_a[ns], offset=off, pos='2', col=strcolors[ns])
      }
    }
}

# Adds a marker to the 2d plot
do_mark = function (marker, marker_error, n, i, j) {

  marker_a = as.numeric(unlist(strsplit(marker, "_")))
  len      = length(marker_a)
  if (len == 0) return
  if (len != n) {
    stop("Wrong length ", len, " of parameters for marker, has to be", n)
  }

  points(marker_a[i], marker_a[j]) 

  marker_error_a = as.numeric(unlist(strsplit(marker_error, "_")))
  len_error = length(marker_error_a)
  if (len_error != n) {
    stop("Wrong length ", len_error, " of parameters for marker_error, has to be ", n)
  }
  if (len != len_error) {
    stop("marker and marker_error: Inconsistent lengths of parameters (", len, len_error, ")")
  }

  x = c(marker_a[i] - marker_error_a[i], marker_a[i] + marker_error_a[i])
  y = c(marker_a[j], marker_a[j])
  lines(x, y)
  x = c(marker_a[i], marker_a[i])
  y = c(marker_a[j] - marker_error_a[j], marker_a[j] + marker_error_a[j])
  lines(x, y)
 
}


# Read command arguments and all existing files
args = commandArgs(TRUE)


log<-file("log_plot_confidence.R")
cat(c("Rscript plot_confidence.R ", args, "\n"), file=log)
close(log)

samples = list()
myplot  = list()
i = 1
ok = T
# Reading all input sample files
while (ok==T) {
  name = args[i]
  if (file.exists(name)) {
    psample = read.table(name)
    samples[[length(samples)+1]] = psample
    ok = T
  } else {
    ok = F
  }
  i = i+1
}

#shade = 1
if (ok==F) {
  if (pmatch("noshade", args[i-1], nomatch=0)!=0) {
    shade = 0
  }
}

nsamples = length(samples)

if (nsamples == 0) {
  stop("Length of sample = 0. Maybe not a valid sample file given on command line?")
}

cat(paste("Plotting", nsamples, "file(s), "))
cat(paste("2d plots: shade = ", shade))
cat(paste(", use_weights = ", use_weights, "\n"))


len = length(fsmooth) 
len1 = len + 1
if (len < nsamples) {
  for (i in len1 : nsamples) {
    fsmooth = append(fsmooth, fsmooth[1])
  }
}


# Getting minimum number of parameters
npar = 10000
for (ns in 1:nsamples) {
  if (npar > dim(samples[[ns]])[2] - 2) {
    npar = dim(samples[[ns]])[2] - 2
  }
}

# ??? rotate y labels
par(las=3)

if (file_test( "-f", configname) == F) {
  stop(c("Configuration file ", configname, " not found"))
}
output = system(paste("get_spar.pl -c ", configname, " R", sep=""), intern=T)
lab    = unlist(strsplit(output, "&"))



# Get parameters from config file
config <- readLines(configname)
for (ln in 1:length(config)) {

  if (pmatch("min", config[ln], nomatch=0)!=0) {
    tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
    #min = as.double(tmp[2:length(tmp)])
    min = as.double(tmp[2:(npar+1)])
  }

  if (pmatch("max", config[ln], nomatch=0)!=0) {
    tmp = unlist(strsplit(config[ln], "[[:blank:]]+", perl=T))
    max = as.double(tmp[2:(npar+1)])
  }

}


# If plot limits are not given, copy those from config file parameters
if (pmins != "") {
  pmin = as.numeric(unlist(strsplit(pmins, "_")))
  if (length(pmin) != length(min)) {
    stop("Wrong length ", length(pmin), " of parameters pmin, has to be ", length(min))
  }
} else {
  pmin = min
}
if (pmaxs != "") {
  pmax = as.numeric(unlist(strsplit(pmaxs, "_")))
  if (length(pmax) != length(max)) {
    stop("Wrong length ", length(pmax), " of parameters pmin, has to be ", length(max))
  }
} else {
  pmax = max
}




fontf = 1.5 # 1.75		# Font scaling factor
axisf = 1.2		# Axis annotation scaling
#lwd   = 2     		# Line width
#fontf  = 2.5
#axisf  = 2.5

# 1D confidence levels
for (i in 1:npar) {

  if (index_i != -1 && index_i != i-1) next
  
  str = sprintf("%d ", i-1)
  cat(str)

  Min = pmin[i]
  Max = pmax[i]

  if (output_format == "ps") {
  	postscript(tmpname, encoding = "TeXtext.enc", paper="special", width=5, height=5, horizontal=F)
  } else if (output_format == "pdf") {
  	pdf(tmpname, paper="special", width=5, height=5)
  }

  first = T
  for (ns in 1:nsamples) {
    psample = samples[[ns]]

    if (min(psample[,i+2])==max(psample[,i+2])) {
      str = sprintf("[no sample %d] ", ns);   
      cat(str);
      next;
    }

    # New (02/07/2018), use weights
    if (use_weights == TRUE) {
        weights = get_weights(psample)
    } else {
        weights = NULL
    }

    t = density(psample[,i+2], weights=weights)
    if (is.na(Min)) {
      # Deduced parameters
      Min = min(t$x)
      Max = max(t$x)
    }

    if (first==T) { myplot = plot }
    else { myplot = points }

    if (solid == TRUE) { lty = 1 }
    else { lty = ns }

    myplot(t$x, t$y/max(t$y), type="l", xlab=parse(text=lab[i]), ylab="posterior",
           xlim=c(Min,Max), ylim=c(0,1), cex.lab=fontf, cex.axis=axisf, lwd=lwd, lty=lty)


    first = F
  }

  tmp = paste("like1d", i-1, sep="_")
  outname = paste(tmp, output_format, sep=".")
  if (output_format == "ps") {
     cmd = sprintf("add_comment_to_ps.pl -o %s %s", outname, tmpname)
  }
  cmd = paste("mv", tmpname, outname, sep=" ")
  system(cmd)

  
  dev.off()

}

if (npar > 1) {
# 2D confidence levels
for  (i in 1:(npar-1)) {

     xMin = min[i]
     xMax = max[i]
     xLab = lab[i]

  
  for (j in (i+1):npar) {

    if (index_i != -1 && index_i != i-1) next
    if (index_j != -1 && index_j != j-1) next

    str = sprintf("%d%d ", i-1, j-1)
    cat(str)


    yMin = min[j]
    yMax = max[j]
    yLab = lab[j]

    if (output_format == "ps") {
      postscript(tmpname, encoding = "TeXtext.enc", paper="special", width=5, height=5, horizontal=F)
    } else if (output_format == "pdf") {
      pdf(tmpname, paper="special", width=5, height=5)
    }

    # Axes limits

    if (is.na(xMin)) {
      xMin = min(z[,1])
      xMax = max(z[,1])
    }
    if (is.na(yMin)) {
      yMin = min(z[,2])
      yMax = max(z[,2])
    }
    
    first = T

    # Plot shaded image
#    shade = 0
    if (shade==T) {

      cat()

      for (ns in 1:nsamples) {

        z = samples[[ns]][,c(i+2,j+2)]   
          xSmooth = (xMax-xMin)/fsmooth[ns]
          ySmooth = (yMax-yMin)/fsmooth[ns]

        # New (02/07/2018), use weights
        if (use_weights == TRUE) {
            weights = get_weights(samples[ns])
        } else {
            weights = list()
        }

        if (min(z[,1])<max(z[,1]) && min(z[,2])<max(z[,2])) {
          NCOL = ns

          try(HPDregionplotDW(mcmc(z), n=N, h=c(xSmooth, ySmooth), prob=prob, ncol=NCOL, weights=weights,
                              xlab=parse(text=xLab), ylab=parse(text=yLab), Xlim=c(xMin, xMax),
                              Ylim=c(yMin, yMax), pXlim=c(pmin[i], pmax[i]), pYlim=c(pmin[j], pmax[j]),
			      cex.lab=fontf, cex.axis=axisf, lwd=lwd, first=first, shade=T), silent=F)
          first = F

        }

      }
    }

    # Draw contours
    for (ns in 1:nsamples) {
      z       = samples[[ns]][,c(i+2,j+2)]
      if (min(z[,1])<max(z[,1]) && min(z[,2])<max(z[,2])) {
        NCOL = ns
        if (solid == TRUE) { lty = 1 }
        else { lty = ns }

        xSmooth = (xMax-xMin)/fsmooth[ns]
        ySmooth = (yMax-yMin)/fsmooth[ns]

        try( HPDregionplotDW(mcmc(z), n=N, h=c(xSmooth, ySmooth), prob=prob, ncol=NCOL, weights=weights,
                            xlab=parse(text=xLab), ylab=parse(text=yLab), Xlim=c(xMin, xMax),
                            Ylim=c(yMin, yMax), pXlim=c(pmin[i], pmax[i]), pYlim=c(pmin[j], pmax[j]),
                            cex.axis=axisf, lty=lty, lwd=lwd, first=first, shade=F), silent=F)
        first = F
      } else {
        str = sprintf("[no sample %d] ", ns);   
        cat(str);
        next;
      }
    }

    # Keys
    if (with_keys == TRUE) {
      do_legend(keystring, pmin[i], pmax[i], pmin[j], pmax[j], nsamples, no_key_line)
    }

    # Title
   
    title(main = title)
    # Marker
    if (nchar(marker) > 0) {
       if (marker_error == "") {
          # Set default error bar
	  for (k in 1:npar) {
             marker_error = paste(marker_error, (pmax[k]-pmin[k])/20, "_", sep="")
	  }
       }
       do_mark(marker, marker_error, npar, i, j)
    }
    

    dev.off()

    tmp = paste("cont2d", i-1, j-1, sep="_")
    outname = paste(tmp, output_format, sep=".")
    #cmd = sprintf("add_comment_to_ps.pl -o %s %s", outname, tmpname)
    cmd = paste("mv", tmpname, outname, sep=" ")
    system(cmd)

  }
}
}

cat("\n")

#warnings()

# New: Create triangle plot and clean up (from plot_confidence.sh)
if (index_i == -1 && index_j == -1) {
  cmd = paste("all_vs_all.pl -b cont2d -e ", output_format, " -l like1d > all_cont2d.tex", sep="")
  system(cmd)

  if (output_format == "ps") {
    system(paste("ldp.sh all_cont2d.tex -q", sep=""))
  } else if (output_format == "pdf") {
    system("pdflatex all_cont2d.tex")
  }
  system("rm -f all_cont2d.{log,dvi,aux} tmptmp.ps Rplots.pdf")
}

q()

# end
