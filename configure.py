#!/usr/bin/env python

# configure.py
# Martin Kilbinger, Karim Benabed 2009
# Configure file for cosmo_pmc


import os
import getopt
import sys
import re


def replace(key, string, data):
  if key is not None:
    res = re.sub("#+\s*" + string + "\s*(?P<sign>\+?)=.*", string + " \g<sign>= " + str(key), data)
    return res
  else:
    return data


def user_choice(name, default, exist):
  # User input
  keep = enter = abort = -1
  ians = 1
  if exist==True:
    #print "[{num}] Use '{path}' as path to {name}? [default]".format(num=ians,
                                                       #name=name, path=default)
    print "[%d] Use '%s' as path to %s?" % (ians, name, default)
    keep  = ians
    ians += 1
    print "[%d] Enter different path" % ians
  else:
    print "[%d] Enter path" % ians

  #print "[{num}] Enter different path".format(num=ians)
  enter = ians
  ians += 1
  #print "[{num}] Abort".format(num=ians)
  print "[%d] Abort" % ians
  abort = ians
  ians += 1
  print "Answer: [1] ",
  ans = sys.stdin.readline()
  ans = ans.strip()
  if ans is "":
    ans = "1"

  ans = int(ans)
  if ans is keep:
    res = os.path.abspath(default)
  elif ans is enter:
    print "Enter path: ",
    res = sys.stdin.readline()
    res = res.strip()
  elif ans is abort:
    print "Aborting"
    res = -1
  else:
    print "Error: Please type '1', '2' or '3'"
    res = -1

  return res


def get_camb(camb):
  # Look for camb directory if not given as command argument
  if camb is None:
    print "Option --camb not given:"
    exist = False
    for pref in "..", os.environ.get("HOME"), "./CMB":
      if exist==True:
        break
      for cname in "camb", "CAMB", "cosmomc/camb":
        path = pref + "/" + cname
        if os.path.isdir(path):
          exist = True
          break

    camb = user_choice("camb", path, exist)

  return camb


def get_wmap(wmap):
  # Look for wmap likelihood directory if not given as command argument
  if wmap is None:
    print "Option --wmap not given:"
    exist = False
    for pref in "..", os.environ.get("HOME"), "./CMB":
      if exist==True:
        break
      for cname in "WMAP7", "WMAP7":
        path = pref + "/" + cname
        if os.path.isdir(path):
          exist = True
          print "path found: " + path
          break

    wmap = user_choice("wmap", path, exist)

  return wmap


def get_pmclib(pmclib):

  # Look for pmclib if not given as command argument
  if pmclib is None:
    print "Option --pmclib not given:"
    exist = False
    # Look for pmclib_v?.? in directories 'pref'
    for pref in "..", ".", os.environ.get("HOME"):
      if exist==True:
        break
      version_PMCLIB_max = 20
      v = version_PMCLIB_max
      while (v>=0):
	v10 = v/10.0
        path = pref + "/pmclib_v%.1f" % v10
        #print "Testing " + path
        if os.path.isdir(path):
	  #print "Found pmclib in " + path
	  exist = True
          break
        v -= 1

    pmclib = user_choice("pmclib", path, exist)

  return pmclib


def correct_path(path):
  if path == '':
     raise IndexError('path not valid')
  if path[0] is '.':
    return os.path.abspath(path)
  else:
    return os.path.expanduser(path)


def usage(ex):
  print "Usage: configure.py [OPTIONS]"
  print "OPTIONS:"
  print "   --cc PROG          Use PROG as C-compiler (default: gcc)"
  print "   --mpicc PROG       Use PROG as MPI C-compiler (default: mpicc)"
  print "   --ld PROG          Use PROG as linker (default: ld)"
  print "   --ar PROG          Use PORG as archive tool (default: ar)"
  print "   --f90 PROG         Use PROG as fortran compiler (default: ifort). Only needed for camb/WMAP support."
  print
  print "   --pmclib PATH      PATH for pmclib"
  print "   --nicaea PATH      Path for nicaea"
  print
  print "   --gsl_prefix PATH  PATH for gsl include files ('PATH/include/gsl') and libraries ('PATH/lib')"
  print "   --fftw3_prefix PATH"
  print "                      PATH for fftw3 include files ('PATH/include') and libraries ('PATH/lib')"
  print
  print "   --inc_mpi FLAGS    Use FLAGS for mpi compiler (e.g.: '-I/path/to/mpi.h'). To add more than one flag,"
  print "                       enclose white-space separated words with quotes"
  print "   --ldirs_mpi FLAGS  Use FLAGS for mpi linking (e.g.: '-L/path/to/libmpi.a')"
  print
  print "   --cflags CFLAGS    Add CFLAGS to compiler flags"
  print "   --lflags LFLAGS    Add LFLAGS to linker flags"
  print
  print "   -c, --cmb          Add CMB (camb+WMAP) support"
  print "   --camb PATH        Path to camb package (if option '-c' given)"
  print "   --wmap PATH        Path to wmap likelihood package (if option '-c' given)"
  print
  print "   --topo PATH        PATH for topolike library, default: not used"
  print
  print "   --lapackdir PATH   PATH for Intel Math Kernel lapack libraries"
  print "   --guidedir PATH    PATH for Intel Math Kernel guide library"
  print "   --ifcoredir PATH   PATH for Intel Math Kernel ifcore library"
  print "   --cfitsiodir PATH  PATH for cfitsio library"
  print
  print "   --installdir PATH  Base path to copy binaries and executables"
  print
  print "   --debug FLAGS      Add debug flags FLAGS"
  print "   -n, --nosave       Do not save file previous setting in 'Makefile.host.save[n]'"

  if ex>0:
    sys.exit(ex)


def main(argv):

  # Command line arguments
  try:
    opts, args = getopt.getopt(argv[1:], "cnh", ["help", "cc=", "mpicc=", "ld=", "ar=", "f90=",
                                                 "pmclib=", "nicaea=",
                                                 "gsl_prefix=", "fftw3_prefix=", "fftw_prefix=",
                                                 "inc_mpi=", "ldirs_mpi=",
                                                 "cflags=", "lflags=",
                                                 "lapackdir=", "guidedir=", "ifcoredir=", "cfitsiodir=",
                                                 "installdir=",
                                                 "cmb", "camb=", "wmap=", "topo=",
                                                 "debug=", 
                                                 "no_save"])
  except getopt.GetoptError, err:
    print err
    usage(1)

  cc = mpicc = ld = ar = f90 = gsl = fftw = None
  pmclib = nicaea = None
  inc_mpi = ldirs_mpi = None
  cflags = lflags = None
  lapackdir = guidedir = ifcoredir = cfitsiodir = installdir = None
  debug = None
  cmb  = 0
  camb = wmap = topo = None
  nosave = 0

  for opt, arg in opts:
    if opt == "--cc":
      cc = arg
    elif opt == "--f90":
      f90 = arg
    elif opt == "--mpicc":
      mpicc = arg
    elif opt == "--ld":
      ld = arg
    elif opt == "--ar":
      ar = arg

    elif opt == "--pmclib":
      pmclib = arg
    elif opt == "--nicaea":
      nicaea = arg

    elif opt == "--gsl_prefix":
      gsl = os.path.abspath(arg)
    elif opt in  ("--fftw_prefix", "--fftw3_prefix"):
      fftw = os.path.abspath(arg)

    elif opt == "--inc_mpi":
      inc_mpi = arg
    elif opt == "--ldirs_mpi":
      ldirs_mpi = arg

    elif opt == "--cflags":
      cflags = arg
    elif opt == "--lflags":
      lflags = arg

    elif opt == "--lapackdir":
      lapackdir = arg
    elif opt == "--guidedir":
      guidedir = arg
    elif opt == "--ifcoredir":
      ifcoredir = arg
    elif opt == "--cfitsiodir":
      cfitsiodir = arg
    elif opt == "--installdir":
      installdir = arg

    elif opt in ("-c", "--cmb"):
      cmb = 1
    elif opt == "--camb":
      camb = arg
    elif opt == "--wmap":
      wmap = arg
    elif opt == '--topo':
      topo = arg
    elif opt == "--debug":
      debug = arg
    elif opt in ("-n", "--no_save"):
      nosave = 1
    elif opt in ("-h", "--help"):
      usage(2)                  
    else:
      print "Error: unhandled option '" + opt + "'"


  mkhost = "Makefile.host"

  saved = 0
  if os.path.exists(mkhost) and nosave is not 1:
    # Save previous file
    count = 1
    while True:
      #print str(count)
      save = "Makefile.host.save%d" % count
      if not os.path.exists(save):
        break
      count += 1

    print "File '" + mkhost + "' with previous settings found, saving under '" + save + "'"
    os.system("cp %s %s" % (mkhost, save))
    saved = 1


  # Replace strings in host-specific Makefile
  o = open(mkhost, "w")
  data = open("Makefile.no_host").read()

  # Write command line to Makefile
  str = ""
  for arg in sys.argv:
    str += arg + " "

  str = "# This Makefile was created with the following command:\n# " + str + "\n"
  data = re.sub("#+\s*Makefile.no_host", "# Makefile.host\n" + str, data)
  data = re.sub("#+\s*Use the script.*", "", data)
  data = re.sub("#+.*host-specific flags.", "", data)

  data = replace(cc, "CC", data)
  data = replace(mpicc, "MPICC", data)
  data = replace(ld, "LD", data)
  data = replace(ar, "AR", data)
  data = replace(f90, "F90", data)

  pmclib = get_pmclib(pmclib)
  if pmclib is -1:
    return 1
  path = correct_path(pmclib)
  data = replace(path, "PMCLIB", data)

  data = replace(nicaea, "NICAEA", data)
 
  data = replace(cmb, "DOWMAP", data)
  if cmb is 0:
    if camb is not None:
      print "Option '--camb' doesn't make sense without '-c' (CMB support)"
      return 1
    if wmap is not None:
      print "Option '--wmap' doesn't make sense without '-c' (CMB support)"
      return 1
  else:
    camb = get_camb(camb)
    if camb is -1:
      return 1
    path = correct_path(camb)
    data = replace(path, "MY_CAMB", data)

    wmap = get_wmap(wmap)
    if wmap is -1:
      return 1
    path = correct_path(wmap)
    data = replace(path, "MY_WMAP", data)

  if topo is not None:
    path = correct_path(topo)
    data = replace(path, "TOPO", data)

  data = replace(gsl, "GSL", data)
  data = replace(fftw, "FFTW", data)

  data = replace(inc_mpi, "INC_MPI", data)
  data = replace(ldirs_mpi, "LDIRS_MPI", data)

  data = replace(cflags, "MKCFLAGS", data)
  data = replace(lflags, "MKLFLAGS", data)

  data = replace(lapackdir, "LAPACKDIR", data)
  data = replace(guidedir, "GUIDEDIR", data)
  data = replace(ifcoredir, "IFCOREDIR", data)
  data = replace(cfitsiodir, "CFITSIODIR", data)
  data = replace(installdir, "INSTALLDIR", data)

  data = replace(debug, "DEBUG", data)

  o.write(data)
  o.close()

  print "File '" + mkhost + "' created. "
  print

  print "Now run "
  print
  print "  make clean; make"
  print
  print "If compilation fails, re-run this script with modified command line options"
  print "(see './configure.py -h' for help), or edit '" + mkhost + "' by hand."

   # Check and set env variables
  cosmopmc = os.environ.get("COSMOPMC")
  if not cosmopmc:
    pwd = os.environ.get("PWD")
    print
    print "Optional: Before running the code, for convenience, you can set the environment variable '$COSMOPMC'",
    if not pwd:
      print " to the directory where this script is located"
    else:
      print "."
      print "In csh,tcsh:\n\n  setenv COSMOPMC " + pwd + "\n"
      print "In bash:\n\n  export COSMOPMC=" + pwd + "\n"

  return 0


if __name__=="__main__":
  sys.exit(main(sys.argv))
