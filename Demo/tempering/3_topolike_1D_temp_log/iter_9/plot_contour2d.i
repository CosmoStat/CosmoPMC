include, "likeli.i"
include, "stuff.i"

write, format="%s\n", "2d plots"
write, format="%s\n", "1d plots"
pmaxall = array(double, 1+0, 0+1)
window, display="", hcp="likeli1d_0"
imin = read_chin("./chi2_0", 1, prior=0, quiet=1)
meansig = array(double, 3)
chi2 = chinall(2,)
plot_chi2_1d, chi2, 1, meansig, pmax, mcmc=1, do_fma=1, do_mean=1, do_sigma=[1, 0, 0], nice=1, outputfile="likeli1d_0.txt", text="alpha_topo", Type=1, color=Blue, width=4, norm=0, key=0, key_str="", height=24
pmaxall(0+1,0+1) = pmax
limits, 0.500000, 6.783200
range, 0, 1
l = limits()
x = l(2)-0.05*(l(2)-l(1))
str = swrite(format="%.3f+%.3f-%.3f", meansig(1), meansig(1)-meansig(2), meansig(3)-meansig(1))
 plt, str, x, l(4)-0.05*(0+1), justify="RT", tosys=1, color=Blue
limits, 0.500000, 6.783200
pltitle_height = 24
xytitles, "alpha_topo", "posterior", [-0.01, 0.025]
pltitle, " "
pdf, "likeli1d_0"


f = open("max_post", "w")
spar = array(string, 1+0)
spar(1) = "alpha_topo"
for (i=1; i<=1+0; i++) {
   write, f, format="%2d %15s", i-1, spar(i)
   for (d=1; d<=0+1; d++) {
     write, f, format=" % 10.5g", pmaxall(i,d)
   }
   writeNL, f
}
close, f
quit
