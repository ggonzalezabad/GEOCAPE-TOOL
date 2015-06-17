# ##################################################
# General notes
# Makefile for the SAO Radiative Transfer Tool
# Based on the original GEOCAPE tool by R.J.D. Spurr
# 
# G. Gonzalez Abad April 2013
####################################################

SHELL= /bin/tcsh

# Name of the compiler and the LINK.f statement should be supplied.

# Ordering of modules follows that required for GNU - the linker
# needs the modules in a certain order (a called module must be
# already compiled and linked before).

################################################
# Path variables (User must change all of these)
################################################

# Paths for too subroutines and tool path
TOOL_PATH = ..
DATA_PATH = $(TOOL_PATH)/datareaders/
NC_PATH   = $(TOOL_PATH)/netcdf/
OBJ       = $(TOOL_PATH)/OBJECTS
WARN_PATH = $(TOOL_PATH)/WARN
EXE_PATH  = .

# GC modules paths 
MOD_PATH  = $(TOOL_PATH)/MODULES
GC_OBJ    = $(TOOL_PATH)/GC_OBJECTS
GC_MOD    = $(TOOL_PATH)/GC_MODULES

# To compile the Vlidort v2p6 modules
VLIDORT_PATH      = $(TOOL_PATH)/vlidort_source
VLIDORT_DEF_PATH  = $(VLIDORT_PATH)/vlidort_def
VLIDORT_MAIN_PATH = $(VLIDORT_PATH)/vlidort_main_FO_1p4
VLID_DEF_PATH     = $(VLIDORT_PATH)/vlidort_def
VLID_MAIN_PATH    = $(VLIDORT_PATH)/vlidort_main_FO_1p4
FO_MAIN_PATH      = $(VLIDORT_PATH)/fo_main_1p4
VSUP_PATH         = $(VLIDORT_PATH)/vsup
VLIDORT_OBJ       = $(TOOL_PATH)/VLIDORT_v2p7_OBJECTS
VLIDORT_MOD       = $(TOOL_PATH)/VLIDORT_v2p7_MODULES

#TOOL_COMPILE = g77 -c -Wimplicit -Wall -ffortran-bounds-check -O2
#VLIDORT_COMPILE = g77 -c -Wimplicit -Wall -ffortran-bounds-check -O2
#LAPACK_COMPILE  = g77 -c -ffortran-bounds-check -O2
#LAPACK_NOPT_COMPILE = g77 -c -ffortran-bounds-check
#F77 = g77

#TOOL_COMPILE = gfortran -c -O3 -mcmodel=medium
#VLIDORT_COMPILE = gfortran -c -O3
#LAPACK_COMPILE  = gfortran -c -O3
#LAPACK_NOPT_COMPILE = gfortran -c -O3
#F77 = gfortran -O3
#LIBS = -L/usr/lib64 -lnetcdff

#TOOL_COMPILE = gfortran -c -g  -fbounds-check -Wall -fbacktrace -mcmodel=medium
#VLIDORT_COMPILE = gfortran -c -g -fbounds-check -Wall -fbacktrace #-mcmodel=medium -fbounds-check -Wall -fbacktrace
#LAPACK_COMPILE  = gfortran -c -g  -fbounds-check -Wall -fbacktrace #-fbounds-check -Wall -fbacktrace
#LAPACK_NOPT_COMPILE = gfortran -c -g  -fbounds-check -Wall -fbacktrace #-fbounds-check -Wall -fbacktrace
#F77 = gfortran -g  -fbounds-check -Wall -fbacktrace -mcmodel=medium #-fbounds-check -Wall -fbacktrace

TOOL_COMPILE = gfortran -c -O3 -mcmodel=medium
VLIDORT_COMPILE = gfortran -c -O3
LAPACK_COMPILE  = gfortran -c -O3
LAPACK_NOPT_COMPILE = gfortran -c -O3
F77 = gfortran -O3 -mcmodel=medium

LIBS = -L/usr/lib64 -lnetcdff

#TOOL_COMPILE        = ifort -c -O3 -zero -mcmodel=medium -shared-intel #-g -traceback -check all -warn all
#VLIDORT_COMPILE     = ifort -c -O3 -zero #-g -traceback -check all -warn all
#LAPACK_COMPILE      = ifort -c -O3 -zero #-g -traceback -check all -warn all
#LAPACK_NOPT_COMPILE = ifort -c -O3
#F77                 = ifort    -O3 -zero -shared-intel #-g -traceback -check all -warn all
#LIBS = -L/data/dumbo/gonzalo/lib/Intel/netcdf/lib -lnetcdf -lirc -lc

#TOOL_COMPILE        = pgf90 -c -O3 #-g -C -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -traceback
#VLIDORT_COMPILE     = pgf90 -c -O3 #-g -C -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -traceback
#LAPACK_COMPILE      = pgf90 -c -O3 #-g -C -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -traceback
#LAPACK_NOPT_COMPILE = pgf90 -c 
#F77                 = pgf90    -O3 #-g -C -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -traceback
#LIBS = -L/home/xliu/OzoneFit/hdfutils/PGI_64bit/netcdf/lib -lnetcdf

#################
# Link definition
#################
LINK.f = $(F77)

###############################################
#  OBJECT MODULES for environment and interface
#         External to main package
###############################################
include ./make.objects.tool

# #####################################################
# Define list of source files for v2p6
# (Note: ordering is important because of dependencies)
# Out sourced to make.vlidort
# #####################################################
include ./make.vlidort

#############
# EXECUTABLES
######GC#######
GCTOOL = geocape_tool

VLIDORT_OBJECTS  = $(F90OBJECTS_FO_GC_SOURCES)                           \
                   $(F90OBJECTS_SOURCES)       $(F90OBJECTS_L_SOURCES)   \
                   $(F90OBJECTS_LPS_SOURCES)   $(F90OBJECTS_LPC_SOURCES) \
                   $(SUP_OBJECTS)

GC_OBJECTS       = $(OBJECTS_GEOCAPE_PREP) $(OBJECTS_NETCDF) $(OBJECTS_MODULE)\
                   $(OBJECTS_GEOCAPE_MAIN)

install  : mk_folders mk_links all
distclean: clean_all rm_folders rm_links
all:       $(VLIDORT_OBJECTS) $(GC_OBJECTS) $(GCTOOL).exe

# --------------------
# The main executables
# --------------------
$(GCTOOL).exe: $(VLIDORT_OBJECTS) $(GC_OBJECTS)
			$(LINK.f) -o $(GCTOOL).exe $(VLIDORT_OBJECTS) $(GC_OBJECTS) $(LIBS)
			mv -fv *genmod* $(WARN_PATH); chmod a+x $(GCTOOL).exe

.PHONY : clean

clean_all:
	rm -vf $(OBJ)/*.o               \
           $(GC_MOD)/*.mod          \
           $(GC_OBJ)/*.o            \
           $(VLIDORT_MOD)/*.mod     \
           $(VLIDORT_OBJ)/*.o       \
           $(WARN_PATH)/*__genmod.* \
           $(EXE_PATH)/*.exe
clean_tool:
	rm -vf $(OBJ)/*.o               \
           $(GC_MOD)/*.mod          \
           $(GC_OBJ)/*.o            \
           $(EXE_PATH)/*.exe
clean_v2p6_vlidort:
	rm -vf $(VLIDORT_MOD)/*.mod \
	       $(VLIDORT_OBJ)/*.o

mk_folders:
	mkdir -p $(OBJ) $(VLIDORT_MOD) $(VLIDORT_OBJ) \
             $(WARN_PATH) $(GC_MOD) $(GC_OBJ) 

mk_links:
	ln -v -s -f /data/tempo1/Shared/HITRAN2012 ../geocape_data/HITRAN

rm_folders:
	rm -rf   $(OBJ) $(VLIDORT_MOD) $(VLIDORT_OBJ) \
             $(WARN_PATH) $(GC_MOD) $(GC_OBJ)

rm_links:
	unlink ../geocape_data/HITRAN

echo:
	echo $(F90OBJECTS_SOURCES)

# #########################
# Pattern rules out sourced
# #########################
include ./make.rules
