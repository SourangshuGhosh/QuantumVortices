#======================================================================
# Makefile for compiling the ns3d, 3D Navier-Stokes hydrodynamics
# #======================================================================

####==================================================
## Compiler options

##-------GENERIC compile
FF = mpif90 
FFLAGS = -O3 -fcheck=all -fno-backtrace -fcheck=bounds

####==============================
### CICADA Unice
#FF = mpiifort
#FFLAGS = -O3

####==============================
### LICALLO OCA
#FF = mpiifort
#FFLAGS = -O3

#FF = mpif90
#FFLAGS = -O3
####==============================
##--------IDRIS--Bluegene---------
## xlf90 compiler
#FF = mpixlf90_r
#FFLAGS = -O3 -g -C -qsigtrap
#FFLAGS = -g -qsigtrap -qoptdebug -qfullpath -C \
# -qflttrap=overflow:underflow:nanq:zerodivide:invalid:qpxstore:enable \
# -qfloat=nans -qsigtrap -qkeepparm
#
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  ## IDRIS--Bluegene/Q Turing
##++++++++++++GCC-Compiler+++++
#FF = mpif90 
#FFLAGS = -O3 -fcheck=all -fno-backtrace -fcheck=bounds

####==============================
### OCCIGEN CINES
#FF = mpiifort
#FFLAGS = -O3 #-check all -CB -traceback #-g
####=====================================================


OBJ1= driver.o
MOD1= FFTW3.o

compile= $(FF) $(FFLAGS) -cpp

####==============================
## FFTW libraries.
####==============================
##-------inca macbookair----------
#link = -lfftw3_mpi -lfftw3 -lm -L/Users/shukla/softwares/fftw/fftw-3.3.8/lib
##------INPHYNI desktop Nice------
#link = -lfftw3_mpi -lfftw3 -lm -L/home/svishwanath/softwares/fftw/fftw-3.3.8/lib
##------IIT KGP desktop------
link = -lfftw3_mpi -lfftw3 -lm -L/home/vishwanath/softwares/fftw/fftw-3.3.8/lib

#======================================================================

sts.exe:$(OBJ1) $(MOD1)
		$(compile) $(OBJ1) $(MOD1) -o sts.exe $(link)

#======================================================================

driver.o: driver.f90 FFTW3.o 
	$(compile) -c -o driver.o \
	driver.f90  FFTW3.o

# Module files
FFTW3.o: FFTW3.f90
	$(compile) -c -o FFTW3.o \
	FFTW3.f90
#======================================================================

clean:
	rm -f *.o *.mod core *.cmp.f90 *.kmo work.* *.d
total_clean:
	rm -f *.o *.mod core *.exe *.out
clean_out:
	rm -f *.out
#======================================================================
