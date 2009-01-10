# Makefile for GNU Make
#
# Note by unknown:
# set your compiler, i used Intel Fortran Compiler
#
# Note by Torfinn Ottesen, October 5. 2007
# ----------------------------------------------
# No joy using
#    GNU Fortran 95 (GCC) 4.1.2 (Ubuntu 4.1.2-0ubuntu4)
#    OS: Ubuntu 7.04
#    At least one of reported errors is fixed in development version
#    of gfortran.   
#
# Work-around, get Intel compiler: ifort
# 1. Obtain Fortran compiler from Intel as non-commercial user.
#    http://www.intel.com/cd/software/products/asmo-na/eng/282048.htm
#    Intel®Fortran Compiler 10.0, Professional and Standard Editions, for Linux*
# 2. Register and note serial number: XXXX-73PJ7T7U
# 3. Download compiler and installation instructions. 
# 4. Unpack downloaded *.tar.gz file and start install:
#    $  tar -xzf *.tar.gz
#    $  cd <created-directory>
#    $  gedit INSTALL.txt    
#       Find procedure: "Perform the following steps to install the compiler."
#    $  ./install.sh  # starts the install
#
# I chose to install as root using sudo since I am using Ubuntu.
# ./install.sh made the decisions as to what to install and where. 
#
# It was a smooth install. Following this I read:
#     Please refer to /opt/intel/fc/10.0.023/doc/Getting_Started.htm
#     to view instructions on how to invoke the product.
#
# A brief 'ifort' summary for convenience:
#     $ source /opt/intel/fc/10.0.xxx/bin/ifortvars.sh  # set-up
#         Does program set-up:
#         PATH, MANPATH, LD_LIBRARY_PATH, INTEL_LICENSE_FILE
#     $ ifort my_source_file.f90
#     $ man ifort  # require man-page installation
#     $ ifort -help
#
# Put the Makefile into directory ...wafo/exec/linux
#       $ cd .../wafo/exec/linux
#       $ make
#
# Lo and behold! It compiles.
# 
#

# set your compiler, i used Intel Fortran Compiler 
F90=ifort

# set your compiler param 
PARAM=-O -fixed

SOURCE_DIR=../../source

TARGET_RIND51=$(addsuffix 51.exe,sp2Acdf sp2mM sp2tccpdf)
TARGET_RIND70=$(addsuffix 70.exe,sp2Acdf sp2mM sp2tccpdf sp2thpdf)

.PHONY: all clean clean-mod clean-all

all: $(TARGET_RIND51) minmax.exe $(TARGET_RIND70) rindd70.exe sp2mmt70.exe simduff.exe down2cc.exe 


%51.exe:	$(SOURCE_DIR)/rind51.f $(SOURCE_DIR)/%1.f 
	$(F90)  $(PARAM) -o $@ $^ 

%70.exe:	$(SOURCE_DIR)/intmodule.f $(SOURCE_DIR)/jacobmod.f $(SOURCE_DIR)/rind70.f $(SOURCE_DIR)/%.f 
	$(F90)  $(PARAM) -o $@ $^ 

sp2mmt70.exe:	$(SOURCE_DIR)/intmodule.f $(SOURCE_DIR)/jacobmod.f $(SOURCE_DIR)/rind70.f $(SOURCE_DIR)/sp2mmtpdf.f 
	$(F90)  $(PARAM) -o $@ $^ 

rindd70.exe:	$(SOURCE_DIR)/intmodule.f $(SOURCE_DIR)/jacobmod.f $(SOURCE_DIR)/rind70.f $(SOURCE_DIR)/rindinterf.f 
	$(F90)  $(PARAM) -o $@ $^ 

minmax.exe:	$(SOURCE_DIR)/rind.f $(SOURCE_DIR)/twog.f $(SOURCE_DIR)/minmax.f
	$(F90)  $(PARAM) -o $@ $^
       
%.exe:	$(SOURCE_DIR)/%.f
	$(F90)  $(PARAM) -o $@ $^ 

clean:
	rm *.exe

clean-mod:
	rm *.mod
	
clean-all:
	rm *.exe *.mod

