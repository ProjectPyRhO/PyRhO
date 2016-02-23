#!/bin/bash
set -e # exit on any error

##### NEURON installation script for Debian Linux (Ubuntu) / OSX #####

# TODO
# Add exports to profile
# Ascertain correct version of xorg-server

useRepo=true
freshTar=true
Ndir=${1:-"$HOME/NEURON"} # e.g. "/home/username/NEURON" 
Ndir=$(echo $Ndir | sed 's:/*$::') # Remove trailing forward-slashes
NRNver=${2:-"nrn-7.4"}
IVver=${3:-"iv-19"}
NRNv="v7.4"
IVdir="iv"
NRNdir="nrn"
pyVer="python3" #"dynamic" # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=6&t=3386
#pyVer=${2:-"3"} # Allow choice of Python version when backported

if [[ -L $0 ]] ; then
	setupDir=$(dirname $(readlink -f $0)) ;
else
	setupDir=$(dirname $0) ;
fi ;

# Run the script without sudo or references to $USER become root
# ./setup.sh /home/user/NEURON
echo "Running $0 from $setupDir..."

if [[ ! -d "${Ndir}" ]]; then
	mkdir ${Ndir}
	#cp $setupDir/mk_hocusr_h.py ${Ndir}
fi


# Linux (Ubuntu) notes
# ====================
# The package installs to /usr/local/bin/
# To remove the package: sudo apt-get remove nrn-7.3.x86_64

# Mac OS X notes
# ==============
# http://www.neuron.yale.edu/neuron/download/compilestd_osx
# See: https://trac.macports.org/attachment/ticket/40457/Portfile-neuron
# http://compneuro.umb.no/LFPy/information.html

# General installation notes
# http://www.davison.webfactional.com/notes/installation-neuron-python/
# http://journal.frontiersin.org/Journal/10.3389/neuro.11.001.2009/full#h9

OS=$(uname)		# Get OS type. 				Alternative: echo $OSTYPE
CPU=$(uname -m)	# Get machine architecture.	Alternative: echo $HOSTTYPE
echo "Operating System ${OS} detected on an ${CPU} architecture!"

echo " *** Installing dependencies... *** "
case $OS in # if [[ "OS" == "Linux" ]]; then
'Linux') # http://stackoverflow.com/questions/394230/detect-the-os-from-a-bash-script
	# http://www.neuron.yale.edu/neuron/download/compile_linux
	# Have more conditionals for Debian vs Red Hat etc. 
	xPath="/usr"
	if [[ -f /etc/debian_version ]]; then
		sudo apt-get update
		# Use \ continuation character for fewer calls to apt-get
		# Consider -y --force-yes for automated installations
		sudo apt-get -y install autotools-dev autoconf automake libtool bison flex
		# bison and flex (and autoconf automake libtool?) are needed if useRepo
		# build-essential python-dev
		
		sudo apt-get -y install xfonts-100dpi libncurses5-dev libxext-dev libreadline-dev
		
		sudo apt-get -y install libopenmpi-dev openmpi-bin openmpi-doc openmpi-common
		#sudo apt-get install python3-mpi4py #python-mpi4py # openmpipython

		#sudo apt-get install git git-man git-doc git-gui gitweb
		#sudo apt-get install python-setuptools python-setuptools-doc python-virtualenv python-pip 
		
		# Fortran compiler for scientific stack i.e. building scipy with pip3. Alternatively just install scipy, numpy and matplotlib
		#sudo apt-get remove g77 # Remove g77 to avoid conflicts with gfortran
		sudo apt-get -y install gcc gfortran liblapack-dev libblas-dev libatlas-base-dev libatlas-dev
		# See here for more details on linear algebra libraries: http://docs.scipy.org/doc/numpy-1.10.1/user/install.html
		# http://docs.scipy.org/doc/scipy/reference/tutorial/linalg.html
		# Alternative Fortran compiler: g77. Do not mix compilers! Check compiler with ldd /usr/lib/{libblas.so.3,liblapack.so.3}
		# TODO: Add option to compile LA library from source for optimal performance
		# export BLAS=/usr/lib/libblas.so
		# export LAPACK=/usr/lib/liblapack.so.3
		# export ATLAS=/usr/lib/libatlas.so
		
		sudo apt-get -y install libfreetype6-dev libxft-dev # Fix a bug in upgrading matplotlib with pip3
		#sudo apt-get build-dep matplotlib
		#sudo pip3 install matplotlib
		
		# NEURON requires Python 2.7 (3 will not work)
		sudo apt-get -y install python3-setuptools python3-pip
		sudo apt-get -y install python3-numpy python3-scipy python3-matplotlib
		sudo pip3 install -U pip
		
		# Install IPython
		# http://ipython.org/ipython-doc/2/install/install.html
		# disutils2 replaces setuptools and pip replaces easy_install
		#pip install ipython[all]
		#pip3 install jupyter # Installed in pip setup.py
		
		if [[ "$useRepo" == true ]] ; then
			sudo apt-get -y install mercurial mercurial-common
		fi
		
		# Install pandoc - a dependency of nbconvert
		# sudo apt-get install pandoc
	fi
	;;
    
'Darwin') # elif [[ "OS" == "Darwin" ]]; then
	# http://www.neuron.yale.edu/neuron/download/compilestd_osx
	xPath="/opt/local"	
	# Assumes that Xcode, Python and IPython are already installed. 
	#sudo port selfupdate
	#sudo port upgrade outdated
	sudo port install automake autoconf libtool xorg-libXext ncurses mercurial bison flex readline
	sudo port install openmpi
	
	### Setup paths to required libraries and tools ###
	# export CC=/opt/local/bin/gcc-mp-4.9
	# export CXX=/opt/local/bin/g++-mp-4.9
	# Alternatively use gcc_select or below
	#sudo port install gcc49
	#sudo port install mpich-gcc49
	#sudo port select --set gcc mp-gcc49
	sudo port select --set mpi openmpi-mp-fortran
	export ACLOCAL_FLAGS="-I/opt/local/share/aclocal $ACLOCAL_FLAGS" # Autoconfigure tools
	export MPICC=/opt/local/bin/mpicc-openmpi-mp #openmpicc
	echo "export MPICC=/opt/local/bin/mpicc-openmpi-mp" >> ~/.profile
	export MPICXX=/opt/local/bin/mpicxx-openmpi-mp #openmpic++
	echo "export MPICXX=/opt/local/bin/mpicxx-openmpi-mp" >> ~/.profile
	export MPIHOME=/opt/local
	echo "export MPIHOME=/opt/local" >> ~/.profile
	#sudo port install xorg-server
	sudo port install xorg-server-devel
	# Determine which X11 version with "defaults read | grep X11"
	#defaults write org.macosforge.xquartz.X11 wm_ffm -bool true # To click on a button without first having to make the window active
	#defaults write org.x.X11 wm_ffm -bool true
	defaults write org.macports.X11 wm_ffm -bool true
	#defaults write org.macports.X11 wm_click_through -bool true # Disables the default behaviour of swallowing window-activated mouse events
	#export PATH=/usr/X11:/usr/X11R6:$PATH
	export LD_LIBRARY_PATH=/opt/local/lib:/opt/local/libexec:$LD_LIBRARY_PATH # Path the linker searches for linking dynamic shared libraries
	echo "export LD_LIBRARY_PATH=/opt/local/lib:/opt/local/libexec:$LD_LIBRARY_PATH" >> ~/.profile
	export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
	echo "export DYLD_LIBRARY_PATH=/opt/local/lib:$DYLD_LIBRARY_PATH" >> ~/.profile
	#export LIBTOOLIZE=/opt/local/bin/glibtoolize
	;;
	
# cygwin? # http://www.neuron.yale.edu/neuron/download/compile_mswin
*) # else
	echo "Unknown OS"
	exit 1
	;;
esac # fi

echo -e " *** Dependencies installed! *** \n"


echo -e "\n************************************************************"
echo "Installing NEURON to $Ndir..."
echo -e "************************************************************\n"

mkdir -p $Ndir
cd $Ndir

if [[ "${freshTar}" == true ]] ; then
	if [[ -d "${IVdir}" ]]; then
		rm -R ${IVdir}
	fi	
	if [[ -d "${NRNdir}" ]]; then
		rm -R ${NRNdir}
	fi
fi


echo " *** Checking for tar files... *** "
case $OS in 
'Linux')
	if [[ "$useRepo" == true ]] ; then
		# N.B. Mercurial will never support Python 3.0 - 3.4!!!
		sudo update-alternatives --install /usr/bin/python python /usr/bin/python2 2 # Ensure python2 is added to update-alternatives before setting it to default		
		sudo update-alternatives --auto python
		hg clone http://www.neuron.yale.edu/hg/neuron/nrn
		hg clone http://www.neuron.yale.edu/hg/neuron/iv
		# Change the prefixes for configure when building from hg repo
	else
		if [[ ! -a $IVver.tar.gz ]]; then
			wget http://www.neuron.yale.edu/ftp/neuron/versions/$NRNv/$IVver.tar.gz
		fi
		if [[ ! -a $NRNver.tar.gz ]]; then
			wget http://www.neuron.yale.edu/ftp/neuron/versions/$NRNv/$NRNver.tar.gz
		fi
	fi
	;;
'Darwin') # Use curl on OSX or install wget with macports
	if [[ "$useRepo" == true ]] ; then
		# N.B. Mercurial will never support Python 3.0 - 3.4!!!
		sudo port select --set python python27
		hg clone http://www.neuron.yale.edu/hg/neuron/nrn
		hg clone http://www.neuron.yale.edu/hg/neuron/iv
		# Change the prefixes for configure when building from hg repo
		sudo port select --set python python34
	else	
		if [[ ! -a $IVver.tar.gz ]]; then
			curl -O http://www.neuron.yale.edu/ftp/neuron/versions/$NRNv/$IVver.tar.gz		
		fi
		if [[ ! -a $NRNver.tar.gz ]]; then
			curl -O http://www.neuron.yale.edu/ftp/neuron/versions/$NRNv/$NRNver.tar.gz
		fi
	fi
	;;
*) # else
	echo "Unknown OS"
	exit 1
	;;
esac
echo -e " *** Required tar files are ready! *** \n"



### Hack!!!
#sudo ln -s /opt/local/bin/glibtoolize /usr/bin/glibtoolize # Needed for build.sh


if [[ "$pyVer" == "dynamic" ]]; then
	if [[ "$OS" == "Linux" ]]; then
		# update-alternatives --list python
		sudo update-alternatives --install /usr/bin/python python /usr/bin/python2 2
		sudo update-alternatives --install /usr/bin/python python /usr/bin/python3 1
		# python --version
		# update-alternatives --list python
		# ls -al /usr/bin/python*
		sudo update-alternatives --set python /usr/bin/python3
	
	elif [[ "$OS" == "Darwin" ]]; then
		sudo port select --set python python34
	fi
fi

# Build and install Interviews
if [[ "${useRepo}" == false ]] ; then #[[ ! -d "${IVdir}" || "${freshTar}" == true ]]; then
	#mkdir $Ndir/$IVdir
	tar xzf $IVver.tar.gz #-C $Ndir/$IVdir # Create dir first. Absolute path?
	if [[ -d "${IVver}" ]]; then	
		mv $IVver $IVdir
	fi
fi

# mv $IVver.tar.gz iv # Rerun with this change carried through
cd $IVdir #cd iv* #
echo -e "\n *** Configuring Interviews... *** "

if [[ "${useRepo}" == true || "$OS" == "Darwin" ]]; then # Is this necessary for OS X when building from tars?
	#aclocal
	#autoheader
	#autoconf
	./build.sh
fi ;

#./configure --prefix=`pwd`
./configure --prefix=$Ndir/$IVdir --with-x --x-includes=$xPath/include/ --x-libraries=$xPath/lib/ 
echo -e " *** Configured Interviews! *** \n"

echo -e "\n *** Making Interviews... *** "
make
echo -e " *** Made Interviews! *** \n"

echo -e "\n *** Installing Interviews... *** "
make install
echo -e " *** Installed Interviews! *** \n"

cd $Ndir
# Build and install NEURON
if [[ "${useRepo}" == false ]] ; then #[[ ! -d "${NRNdir}" || "${freshTar}" == true ]]; then
	#mkdir $Ndir/$NRNdir
	tar xzf $NRNver.tar.gz #-C $Ndir/$NRNdir
	if [[ -d "${NRNver}" ]]; then	
		mv $NRNver $NRNdir
	fi	
fi


# mv $NRNver.tar.gz iv # Rerun with this change carried through
cd $NRNdir #cd nrn* #
echo -e "\n *** Configuring NEURON... *** "

if [[ "$useRepo" == false ]] ; then
	sudo sh ./src/nrnmpi/mkdynam.sh # make dependency bug: http://www.neuron.yale.edu/phpbb/viewtopic.php?f=4&t=3049
fi

if [[ "${useRepo}" == true || "$OS" == "Darwin" ]]; then # Is this necessary for OS X when building from tars?
	./build.sh
fi ;

### Hacks!!! Run these if Python 3 is the system default ###
### Alternatively set the system version to python 2: sudo port select python python27
### then link configure script to python 3: --with-nrnpython=/opt/local/bin/python3.4
echo -e "\n\n ***** Hacking broken files!!! *****\n\n"
cd $Ndir
sudo chown -R $USER . #$Ndir

# Replace the following hacks with: 2to3 -w mk_hocusr.h
# https://www.neuron.yale.edu/phpBB/viewtopic.php?f=6&t=3386&p=14346#p14346
#cp mk_hocusr_h.py $NRNdir/src/oc/. # Changed print commands
2to3 -w $NRNdir/src/oc/mk_hocusr_h.py
sed -i '1i from __future__ import print_function' $NRNdir/src/oc/mk_hocusr_h.py

#cp configure $NRNdir/. # Changed print statement on line 6599
if [[ "$OS" == "Darwin" ]]; then
	sudo sed -i .bak -e "s/print sys.api_version,/from __future__ import print_function; print(sys.api_version)/" $NRNdir/configure
else # different syntax for -i on linux
	#sudo sed -i.bak -e "s/print sys.api_version,/print(sys.api_version)/" $NRNdir/configure
	sudo sed -i.bak -e "s/print sys.api_version,/from __future__ import print_function; print(sys.api_version)/" $NRNdir/configure
fi

cd $NRNdir
#############################################################


# Could change --with-nrnpython=dynamic to --with-nrnpython=python3
#./configure --prefix=$Ndir/$NRNdir --with-iv=$Ndir/$IVdir --with-nrnpython=dynamic --with-paranrn=dynamic \
#--with-x --x-includes=$xPath/include/ --x-libraries=$xPath/lib/
./configure --prefix=$Ndir/$NRNdir --with-iv=$Ndir/$IVdir --with-nrnpython=$pyVer --with-paranrn=dynamic \
--with-x --x-includes=$xPath/include/ --x-libraries=$xPath/lib/
#./configure --prefix=$Ndir/$NRNdir --without-iv --with-nrnpython=dynamic --with-paranrn=dynamic \
#--with-x --x-includes=$xPath/include/ --x-libraries=$xPath/lib/
# --with-nrnpython=/opt/local/bin/python
echo -e " *** Configured NEURON! *** \n"

echo -e "\n *** Making NEURON... *** "
make
echo -e " *** Made NEURON! *** \n"

echo -e "\n *** Installing NEURON... *** "
make install
echo -e " *** Installed NEURON! *** \n"

export PATH=$Ndir/$IVdir/$CPU/bin:$PATH
export PATH=$Ndir/$NRNdir/$CPU/bin:$PATH # Necessary?
export PATH=$Ndir/$CPU/bin:$PATH

echo "export PATH=$Ndir:$Ndir/$NRNdir/$CPU/bin:$Ndir/$IVdir/$CPU/bin:$PATH" >> ~/.profile # :$Ndir/$CPU/bin

#export LD_LIBRARY_PATH=/opt/local/lib/
#source $Ndir/nrnenv

# Build the icons for GUI use
cd $Ndir/$NRNdir

if [[ "$OS" == "Darwin" ]]; then
	make after_install
fi

echo -e "\n *** Setting up nrnpython... ***"
cd src/nrnpython
if [[ -x "$(command -v python3)" ]]; then # Explictly use Python 3
	sudo python3 setup.py install
else # Assume python already points to python3
	sudo python setup.py install # Look at installoption
fi

#export PYTHONPATH=$Ndir/lib/python:$PYTHONPATH
# Create environment variables so python can import the NEURON module
export NEURONHOME=$Ndir/$NRNdir #/NEURON-7.3/nrn
export PYTHONPATH=$PYTHONPATH:$NEURONHOME/lib/python
echo "export PYTHONPATH=$PYTHONPATH:$NEURONHOME/lib/python" >> ~/.profile
export NRN_NMODL_PATH=$Ndir
echo "export NRN_NMODL_PATH=$Ndir" >> ~/.profile
#cd ~
source ~/.profile
echo -e " *** Set up nrnpython! *** \n"

cd $Ndir

### Run tests ###
# idraw
# nrngui
# neurondemo

# Parallel tests
# cd $Ndir/$NRNdir/src/parallel
# mpirun -n 4 $N/$CPU/bin/nrniv -mpi test0.hoc
# mpirun -n 4 $N/$CPU/bin/nrniv -mpi test0.py

### Download and Compile mod files ###
# nrnivmodl
