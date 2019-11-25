#!/usr/bin/env bash

# :Name: install_CosmoPMC.sh
# :Author: Martin Kilbinger
# (adapted from install_shapepipe by Sam Farrens)
# :Date: 2019
# :Package: CosmoPMC
# :Description: Install CosmoPMC and dependend libraries.


##################
# Global variables
##################

# Divider line
line="########################################################################"

# Version
version=0.3

# Default values for variables
BASE_DIR=$PWD
BUILD_DIR=$BASE_DIR/build
PMCENV="cosmopmc"
TOPO=""
LFLAGS=""

# Help string
help="$(basename "$0") [OPTIONS]\n\n
Options:\n
\t-h, --help\t\t show this help message and exit\n
\t--build-dir PATH\t set CosmoPMC build path (default is \$PWD/build)\n
\n
Executable and library options:\n
\t--lflags LFLAGS\t\t linker flags\n
\t--topo PATH\t\t topology likelihood path (default not used)\n
\n
"

#############
# Subroutines
#############

# Print start message
start() {
  echo ''
  echo $line
  echo 'CosmoPMC installation script'
  echo ''
  echo 'Author: Martin Kilbinger'
  echo 'Year: 2019'
  echo 'Version:' $version
  echo $line
  echo ''
}

# Check if a binary executable is already installed in the conda environment
check_conda() {
  if ! type -t "conda" > /dev/null
  then
    echo "Conda command not found, make sure it is installed before proceding."
    exit
  else
    echo "Conda command found"
  fi
} 

print_setup() {
  echo 'Operating system: ' $SYSOS
  echo 'Base directory (CosmoPMC clone): ' $BASE_DIR
  echo 'Build directory: ' $BUILD_DIR
  echo 'Environment name: ' $PMCENV
  echo 'Topolike directory: ' $TOPO
  echo 'Linker flags: ' $LFLAGS
  echo ''
}

# Function to report progress
report_progress() {
  echo ''
  echo $line
  echo 'Installing' $1
  echo $line
  echo ''
}



##############
# Start script
##############

# Parse command line options
for i in "$@"
do  
case $i in
    -h|--help)
    start
    echo -ne $help
    shift
    exit
    ;;
    --build-dir=*)
    BUILD_DIR="${i#*=}"
    shift
    ;;
    --lapackdir=*)
    TOPO="${i#*=}"
    shift
    ;;
    --topo=*)
    TOPO="${i#*=}"
    shift
    ;;
esac
done


## Create conda environment

# Start script
start

# Check if conda is installed
check_conda

# Find the operating system 
case "$OSTYPE" in
  darwin*)
  SYSOS="macOS"
  ;;
  linux*)
  SYSOS="LINUX"
  ;;
  *)
  echo "unknown: $OSTYPE"
  exit 1
  ;;
esac

# Create build directory if it does not already exist
if [ ! -d "$BUILD_DIR" ]
then
  mkdir $BUILD_DIR
fi

# Print script set-up
print_setup

# Build conda environment
conda info --envs | grep $PMCENV > /dev/null
x=$?
if [ $x == 0 ]; then
  echo "CosmoPMC environment already exists, skipping create"
else
  report_progress 'CosmoPMC environment'
  conda env create -f environment.yml
fi

# Activate conda environment
if [ "$CONDA_DEFAULT_ENV" == "$PMCENV" ]; then
  echo "$PMCENV conda environment already active, skipping activate"
else
  #eval "$(conda shell.bash hook)"
  source activate $PMCENV
  echo "Activating conda environment" $CONDA_PREFIX
fi

############################
# install macos requirements
############################

# Set up macOS environment
if [ "$SYSOS" == "macOS" ]; then
  report_progress 'macOS Requirements'
  conda install -n shapepipe wget -y
  conda install -n shapepipe automake autoconf libtool -y
  export C_INCLUDE_PATH=$CONDA_PREFIX/include
  export CFLAGS="-Wl,-rpath,$CONDA_PREFIX/lib"
fi

##########################
# install mpi requirements
##########################

conda install -n $PMCENV -c conda-forge mpich -y

##########################
# Build external libraries
##########################

report_progress 'fftw3'
conda install -n $PMCENV -c conda-forge fftw -y

report_progress 'gsl'
conda install -n $PMCENV -c conda-forge gsl=1.16 -y

report_progress 'lacpack'
conda install -n $PMCENV -c conda-forge liblapack -y


# C-compiler stuff
if [ "$SYSOS" == "LINUX" ]; then
  report_progress 'gxx_linux-64'
  conda install -n $PMCENV gxx_linux-64 -y
fi

#########################
# Build external programs
#########################

report_progress "cmake"
conda install -n cosmopmc -c conda-forge cmake -y

report_progress "gnuplot"
conda install -n cosmopmc -c conda-forge gnuplot -y

#############################
# Build PMC-related libraries
#############################

# pmclib
report_progress 'pmclib'
cd $BUILD_DIR
if [ -e pmclib ]; then
  rm -rf pmclib
fi
git clone https://github.com/CosmoStat/pmclib.git
cd pmclib
python -c 'import sys; print(open("wscript").read().replace("MK_TO_REPLACE", sys.argv[1]))' $CONDA_PREFIX > wscript.new
mv wscript.new wscript
python2 ./waf configure --prefix=$CONDA_PREFIX --m64
python2 ./waf build install
cd $BASE_DIR

# nicaea
report_progress 'nicaea'
cd $BUILD_DIR
if [ -e nicaea ]; then
  rm -rf nicaea
fi
git clone https://github.com/CosmoStat/nicaea.git
mkdir nicaea/build
cd nicaea/build
cmake .. -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
make && make install
cd $BASE_DIR

# CosmoPMC
report_progress 'CosmoPMC'
cd $BASE_DIR

if [ "$TOPO" == "" ]; then
  arg_topo=""
else
  arg_topo="--topo $TOPO"
fi
if [ "$LFLAGS" == "" ]; then
  lflags=""
else
  lflags="$LFLAGS"
fi

python2 ./configure.py --pmclib=$CONDA_PREFIX --nicaea=$CONDA_PREFIX --installdir=$CONDA_PREFIX $arg_topo $LFLAGS  # --inc_mpi -I/usr/include/mpich-x86_64 --ldirs_mpi=-L/usr/lib64/mpich/lib # --lflags -lm
make && make install
