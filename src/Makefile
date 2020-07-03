#======================================================================
# Makefile for compiling the 2D hydrodynamics
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

#FF = h5pfc
#FFLAGS = -O3 -g -traceback -check all -fp-stack-check -CB

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
#FF = h5pfc
#FFLAGS = -O3 #-check all -CB -traceback #-g
####==============================
### SNBOSE cluster PHYS-IITKGP
#FF = h5pfc
#FF = mpifort
#FFLAGS = -O3 #-check all -CB -traceback #-g
####=====================================================

OBJ1= driver.o initial_config.o solve.o pspectral.o \
	   diagnostics.o io_fields.o gpe.o gpe_particle.o
MOD1= mod_precision.o mod_fft.o mod_constants.o mod_2dHD.o mod_grid.o \
	   FFTW3.o f90mtranu.o mod_mpi.o mod_bench.o mod_gpparticles.o

###==================================
## Conditional compilation keys
##-----------------------------

### Forcing options
#IOSCHEME = -DIOHDF5
IOSCHEME = -DIOEVERYPROC
OUTPUTSP = -DFOURSPACE
OUTPUTSP = #-DPHYSSPACE
STEPSCHEME = -DFIXDT
DTVAL = #-DDT_VAL_FROM_FILE
## adaptive time scheme is crap here.
#STEPSCHEME = -DADAPT
#HD = $(FLUID2) $(STEPSCHEME) $(FORCED) #$(SCALAR) #$(FLUID2)
GPADV = #-DGPE_ADVECTION
GPGL = #-DGPE_ARGLE $(GPADV)
GPCOUNTFLOW = #-DGPE_COUNTERFLOW
GPSELFTRUNC = #-DGPE_SELFTRUNC
GPSGL = #-DGPE_SGLE $(GPSELFTRUNC) $(GPCOUNTFLOW)
GPATRCT = #-DGPE_ATTRACT
GPFORCED = #-DGPE_FRCDISP
GPPARTICLE = #-DGPPART
GPSTSPECTRA = #-DGPSTS
GPOTOCS = -DGPOTOC
HD = -DSOL_GPE
CPPFLAGS = $(HD) $(STEPSCHEME) $(GPFORCED) $(GPATRCT) $(GPPARTICLE) $(GPGL) $(GPSGL) $(GPSTSPECTRA) $(IOSCHEME) $(OUTPUTSP) $(GPOTOCS)

####==============================
## Compiler flags collected together
####==============================
##%%%%%%%%%%General purpose%%%%%%%
compile= $(FF) $(FFLAGS) $(CPPFLAGS) -cpp


####==============================
## FFTW libraries.
####==============================
##-------inca macbookair----------
#link = -lfftw3_mpi -lfftw3 -lm -L/Users/shukla/softwares/fftw/fftw-3.3.8/lib

##------PHYSICS IIT KGP------
link = -lfftw3_mpi -lfftw3 -lm -L/home/vishwanath/softwares/fftw/fftw-3.3.8/lib

##-------SNBOSE PHYS IITKGP----------
#link = -lfftw3_mpi -lfftw3 -lm -L/home/vshukla/softwares/fftw/fftw-3.3.8_intel/lib

##------INPHYNI desktop Nice------
#link = -lfftw3_mpi -lfftw3 -lm -L/home/svishwanath/softwares/fftw/fftw-3.3.8/lib

##------numerics ENS paris--------
#link = -lfftw3_mpi -lfftw3 -lm -L/home/vishwanath/softwares/fftw/fftw-3.3.4/lib

##------OLD Desktop CEA-----------
#link = -lfftw3_mpi -lfftw3 -lm -L/home/vshukla/softwares/fftw/fftw-3.3.4/lib

##------CEAWORKSTATION------------
#link = -lfftw3_mpi -lfftw3 -lm -L/opt/fftw/fftw-3.3.7-gfortran/lib

##------MESOPSL PARIS-------------
#link = -lfftw3_mpi -lfftw3 -lm -L/obs/vshukla/softwares/fftw/fftw-3.3.4/lib

##------CICADA Unice--------------
#link = -lfftw3_mpi -lfftw3 -lm -L/softs/fftw/fftw-3.3.3/fftw-3.3.3_intel-12.1/lib

##-------LICALLO------------------
#link = -lfftw3_mpi -lfftw3 -lm -L/beegfs/home/svishwan/softwares/fftw/fftw-3.3.4/lib
#link = -lfftw3_mpi -lfftw3 -lm -L/trinity/shared/OCA/softs/fftw-3.3.7-intel17//lib

##-------Bluegene-----------------
#link = $(WRAPPER_LDFLAGS)
#link = -lfftw3_mpi -lfftw3 -lm -L/bglocal/cn/pub/FFTW/3.3.3/lib
#link = -lfftw3_mpi -lfftw3 -lm -L/linkhome/rech/xhr/rxhr001/fftw/fftw-3.3.4/lib

##-------IDRIS--------------------
#link = -lfftw3_mpi -lfftw3 -lm -L/gpfslocalsup/spack_soft/fftw/3.3.8/intel-19.0.4-ifhmvtovzozeif6roxv5xyn2cmo5grlq/lib

##------Occigen-------------------
#link = -lfftw3_mpi -lfftw3 -lm -L/opt/software/occigen/libraries/fftw/3.3.6-pl2/intel/17.0/intelmpi/2017.0.098/lib
#======================================================================

#======================================================================

vikshobha.exe:$(OBJ1) $(MOD1)
		$(compile) $(OBJ1) $(MOD1) -o vikshobha.exe $(link)

#======================================================================

driver.o: driver.f90 mod_precision.o mod_constants.o mod_2dHD.o mod_fft.o \
	mod_grid.o FFTW3.o  mod_bench.o mod_gpparticles.o
	$(compile) -c -o driver.o \
	driver.f90  mod_precision.o mod_constants.o mod_2dHD.o mod_fft.o \
	mod_grid.o FFTW3.o mod_bench.o mod_gpparticles.o

initial_config.o: initial_config.f90 mod_precision.o mod_constants.o mod_2dHD.o \
	mod_fft.o mod_grid.o FFTW3.o mod_mpi.o f90mtranu.o
	$(compile) -c -o initial_config.o \
	initial_config.f90  mod_precision.o mod_constants.o mod_2dHD.o mod_fft.o \
	mod_grid.o FFTW3.o mod_mpi.o f90mtranu.o

solve.o: solve.f90 mod_precision.o mod_constants.o mod_grid.o f90mtranu.o
	$(compile) -c -o solve.o \
	solve.f90 mod_precision.o mod_constants.o mod_grid.o f90mtranu.o

pspectral.o: pspectral.f90 mod_precision.o mod_constants.o \
	mod_2dHD.o mod_fft.o mod_grid.o FFTW3.o mod_mpi.o
	$(compile) -c -o pspectral.o \
	pspectral.f90 mod_precision.o mod_constants.o \
	mod_2dHD.o mod_fft.o mod_grid.o FFTW3.o mod_mpi.o

diagnostics.o: diagnostics.f90 mod_precision.o mod_constants.o mod_2dHD.o \
	mod_grid.o FFTW3.o mod_fft.o mod_mpi.o mod_gpparticles.o
	$(compile) -c -o diagnostics.o \
	diagnostics.f90 mod_precision.o mod_constants.o mod_2dHD.o \
	mod_grid.o FFTW3.o mod_fft.o mod_mpi.o mod_gpparticles.o

gpe.o: gpe.f90 mod_precision.o mod_constants.o mod_2dHD.o mod_grid.o \
	FFTW3.o mod_fft.o mod_mpi.o f90mtranu.o
	$(compile) -c -o gpe.o \
	gpe.f90 mod_precision.o mod_constants.o mod_2dHD.o mod_grid.o \
	FFTW3.o mod_fft.o mod_mpi.o f90mtranu.o

gpe_particle.o: gpe_particle.f90 mod_precision.o mod_constants.o mod_2dHD.o mod_grid.o \
	FFTW3.o mod_fft.o mod_mpi.o f90mtranu.o mod_gpparticles.o
	$(compile) -c -o gpe_particle.o \
	gpe_particle.f90 mod_precision.o mod_constants.o mod_2dHD.o mod_grid.o \
	FFTW3.o mod_fft.o mod_mpi.o f90mtranu.o	mod_gpparticles.o

io_fields.o: io_fields.f90 mod_precision.o mod_constants.o mod_2dHD.o \
	mod_fft.o mod_grid.o FFTW3.o mod_mpi.o
	$(compile) -c -o io_fields.o \
	io_fields.f90  mod_precision.o mod_constants.o mod_2dHD.o mod_fft.o \
	mod_grid.o FFTW3.o mod_mpi.o

# Module files
mod_precision.o: mod_precision.f90
	$(compile) -c -o mod_precision.o \
	mod_precision.f90

mod_constants.o: mod_constants.f90 mod_precision.o
	$(compile) -c -o mod_constants.o \
	mod_constants.f90 mod_precision.o

mod_grid.o: mod_grid.f90 mod_precision.o FFTW3.o mod_constants.o
	$(compile) -c -o mod_grid.o \
	mod_grid.f90 mod_precision.o FFTW3.o mod_constants.o

mod_fft.o: mod_fft.f90 mod_precision.o mod_constants.o FFTW3.o mod_grid.o mod_2dHD.o mod_mpi.o
	$(compile) -c -o mod_fft.o \
	mod_fft.f90 mod_precision.o mod_constants.o FFTW3.o mod_grid.o mod_2dHD.o mod_mpi.o

mod_2dHD.o: mod_2dHD.f90 mod_precision.o FFTW3.o mod_grid.o
	$(compile) -c -o mod_2dHD.o \
	mod_2dHD.f90 mod_precision.o FFTW3.o mod_grid.o

mod_gpparticles.o: mod_gpparticles.f90 mod_precision.o mod_constants.o mod_grid.o 
	$(compile) -c -o mod_gpparticles.o \
	mod_gpparticles.f90 mod_precision.o mod_constants.o mod_grid.o

FFTW3.o: FFTW3.f90
	$(compile) -c -o FFTW3.o \
	FFTW3.f90

mod_mpi.o: mod_mpi.f90
	$(compile) -c -o mod_mpi.o \
	mod_mpi.f90

mod_bench.o: mod_bench.f90
	$(compile) -c -o mod_bench.o \
	mod_bench.f90

f90mtranu.o: f90mtranu.f90
	$(compile) -c -o f90mtranu.o \
	f90mtranu.f90

#======================================================================

clean:
	rm -f *.o *.mod core *.cmp.f90 *.kmo work.* *.d
total_clean:
	rm -f *.o *.mod core *.exe *.out
clean_out:
	rm -f *.out
#======================================================================
