!!**************************************************************!!
!!                                                              !!
!!                     VIKSHOBHA                                !!
!!                                                              !!
!!--------------------------------------------------------------!!
!! Programming Language: FORTRAN 90                             !!
!! Hydrodynamic solver based on Fourier pseudospectral method.  !!
!! MPI Parallel code with slab decomposition.                   !!
!! Relies on FFTW librarires (Version 3.3.x). 		        !!
!!--------------------------------------------------------------!!
!! Solvers:                                                     !!
!! -------------                                                !!
!! ( I   ) 2D Gross-Pitaevskii equation.                        !!
!!--------------------------------------------------------------!! 
!!                                                              !!
!! Developed by:                                                !!
!!               Vishwanath Shukla                              !!
!!                                                              !!
!! Current affiliation:                                         !!
!!---------------------                                         !!
!! Department of Physics,                                       !!
!! Indian Institute of Technology Kharagpur                     !!
!! Kharagpur, India                                             !!
!! Email 1: vishwanath.shukla@phy.iitkgp.ac.in                  !!
!! Email 2: research.vishwanath@gmail.com                       !!
!!**************************************************************!!

PROGRAM driver

!======================================================================
! Include modules.
!======================================================================

USE FFTW3             !!--Include FFTW3 file. 
USE mod_precision     !!--Sets the precision (single/double) in KIND=GP.
USE mod_grid          !!--Declares and defines space/time resoltuion.
USE mod_bench         !!--For benchmarking.
USE mod_fft           !!--Declares, allocates variables. Sets up FFTs.
USE mod_constants     !!--Defines many constants used in the code.
USE mod_2dHD          !!--Defines important physical parameters + more.
USE mod_mpi
#ifdef GPPART
USE mod_gpparticles
#endif

IMPLICIT NONE


!======================================================================
! List of the input parameters.
!======================================================================
NAMELIST/ckprst/nrst,sprst,strst
NAMELIST/resolution/Nx,Ny
NAMELIST/fileparam/nsplit_orignal,nsplit_current
NAMELIST/tparam/dt_fix,cflmax,maxiter,fsteps,ssteps,gsteps,ststeps
!--------
#ifdef SOL_GPE
NAMELIST/gpparam/xi,csound,rho_avg,disp_factor
NAMELIST/gpfiniteT/oneybeta,nuN,chempot0
#ifdef GPE_COUNTERFLOW
NAMELIST/gpcounterflow/vnx,vny
#endif
#ifdef GPE_SELFTRUNC
NAMELIST/selftruncation/ktrunc
#endif
#ifdef GPE_FRCDISP
NAMELIST/gpforceparam/famp_gp,kf1_gp,kf2_gp
NAMELIST/gphypvis/nhypvis1_gp,nhypvis2_gp,vishypo_gp,vishype_gp,friclowk_gp
#endif
#ifdef GPE_ATTRACT
NAMELIST/gpfocus/gstrength
#endif
#ifdef GPPART
NAMELIST/gpparticle/gp_N_obj,gp_mass_part,ppotd,ppotV0,ppotthick
#endif
#ifdef GPOTOC
NAMELIST/gpotoc/k_otoc_gp,eps_otoc_gp
#endif
#endif
!--------


!======================================================================
!       Setting up the MPI
!======================================================================
!       Initialize
CALL mpi_init(ierr)
! Get the number of processes
CALL mpi_comm_size(mpi_comm_world, nprocs, ierr)
! Get the rank
CALL mpi_comm_rank(mpi_comm_world,myrank, ierr)
! Get the clock time
CALL mpi_barrier(mpi_comm_world,ierr)
time_ini = mpi_wtime()
! Initialize the fftw3
CALL fftw_mpi_init
!======================================================================
!======================================================================


IF (myrank .EQ. 0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!!                  RUN INFORMATION              !!'                 
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!!**************************************************************!!'
PRINT*, '!!                                                              !!'
PRINT*, '!!                     VIKSHOBHA                                !!'
PRINT*, '!!                                                              !!'
PRINT*, '!!--------------------------------------------------------------!!'
PRINT*, '!! Programming Language: FORTRAN 90                             !!'
PRINT*, '!! Hydrodynamic solver based on Fourier pseudospectral method.  !!'
PRINT*, '!! MPI Parallel code with slab decomposition.                   !!'
PRINT*, '!! Relies on FFTW librarires (Version 3.3.x). 	            	 !!'
PRINT*, '!!--------------------------------------------------------------!!'
PRINT*, '!! Solvers:                                                     !!'
PRINT*, '!! -------------                                                !!'
PRINT*, '!! ( I   ) 2D Gross-Pitaevskii equation.                        !!'
PRINT*, '!!--------------------------------------------------------------!!'
PRINT*, '!!                                                              !!'
PRINT*, '!! Developed by:                                                !!'
PRINT*, '!!               Vishwanath Shukla                              !!'
PRINT*, '!!                                                              !!'
PRINT*, '!! Current affiliation:                                         !!'
PRINT*, '!!---------------------                                         !!'
PRINT*, '!! Department of Physics,                                       !!'
PRINT*, '!! Indian Institute of Technology Kharagpur                     !!'
PRINT*, '!! Kharagpur, India                                             !!'
PRINT*, '!! Email 1: vishwanath.shukla@phy.iitkgp.ac.in                  !!'
PRINT*, '!! Email 2: research.vishwanath@gmail.com                       !!'
PRINT*, '!!**************************************************************!!'
PRINT*, '!!==============================================================!!'
PRINT*, '!!                     SOLVER CHOSEN                            !!'
#ifdef SOL_GPE
PRINT*, '!!                   Gross-Pitaevskii Equation                  !!'
#endif
PRINT*, '!!==============================================================!!'
ENDIF


!======================================================================
! Read the input parameters.
!======================================================================

!! Read from the 0th process and then broadcast to all others.
IF (myrank .EQ. 0) THEN

OPEN(UNIT=1,FILE='parameter.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=ckprst,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)

OPEN(UNIT=1,FILE='parameter.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=resolution,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)

OPEN(UNIT=1,FILE='parameter.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=fileparam,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)

OPEN(UNIT=1,FILE='parameter.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=tparam,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)


!--------
#ifdef SOL_GPE
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpparam,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)

OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpfiniteT,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)

#ifdef GPE_COUNTERFLOW
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpcounterflow,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)
#endif

#ifdef GPE_SELFTRUNC
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=selftruncation,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)
#endif

#ifdef GPE_FRCDISP
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpforceparam,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)

OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gphypvis,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)
#endif

#ifdef GPE_ATTRACT
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpfocus,IOSTAT=ioerr,IOMSG=iomesg)
CLOSE(1)
#endif

#ifdef GPPART
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpparticle,IOSTAT=ioerr,IOMSG=iomesg)
#endif

#ifdef GPOTOC
OPEN(UNIT=1,FILE='parameter_gpe.txt',STATUS='UNKNOWN',FORM='FORMATTED')
READ(1,NML=gpotoc,IOSTAT=ioerr,IOMSG=iomesg)
#endif

#endif
!--------

PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! PARAMETERS'
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! MPI ioerr'
PRINT*,ioerr!,iomesg
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! Restart file variables'
PRINT*,nrst,sprst,strst
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! Spatial resolution'
PRINT*, 'Nx:',Nx,'Ny:',Ny
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! Parameters determing data split for writing'
PRINT*, 'Orignal data split into files', nsplit_orignal
PRINT*, 'Processed data split into files', nsplit_current
PRINT*, '!!-----------------------------------------------!!'
PRINT*, 'Time step:', dt_fix
PRINT*, 'Max. of CFL allowed:', cflmax
PRINT*, 'Maxium no. of time steps:', maxiter
PRINT*, 'Fields are stored after fsteps:', fsteps
PRINT*, 'Spectra are stored after ssteps:', ssteps
PRINT*, 'Global data (eg. energy) are stored after gsteps:', gsteps
PRINT*, 'Spatiotemproal spectra (GPE) stored after ststeps:', ststeps
PRINT*, '!!-----------------------------------------------!!'
!!---------------------------------------
#ifdef SOL_GPE
print*, 'GPE params',xi,csound,rho_avg,disp_factor
print*, 'GPE finiteT',oneybeta,nuN,chempot0
#ifdef GPE_COUNTERFLOW
PRINT*, 'CounterFlow vnx and vny', vnx, vny
#endif
#ifdef GPE_SELFTRUNC
PRINT*, 'GPE imposed truncation wave number', ktrunc
#endif
#ifdef GPE_FRCDISP
print*, 'GPE FRC', famp_gp,kf1_gp,kf2_gp
print*, 'GPE DISP',nhypvis1_gp,nhypvis2_gp,vishypo_gp,vishype_gp,friclowk_gp
#endif
!#ifdef GPE_ATTRACT
!print*,
!#endif
#ifdef GPPART
print*, 'GPE Particle', gp_N_obj,gp_mass_part,ppotd,ppotV0,ppotthick
print*, 'Note some parameters of the potential are modified later.'
#endif
#ifdef GPOTOC
PRINT*, 'GPE OTOC, wave number', k_otoc_gp
PRINT*, 'GPE OTOC, Perturbation', eps_otoc_gp
#endif
#endif

!!----------------------------------------------
!! Set the time step for the simulation.       !
!! Note: currently only fixed time step works. ! 
!!----------------------------------------------
dt = dt_fix
dt_adpt_ns = dt_fix
dt_adpt_2f = dt_fix
dt_adpt_scl1 = dt_fix

ENDIF

CALL MPI_BCAST(Nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(Ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nsplit_orignal,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nsplit_current,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(dt_fix,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(cflmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(maxiter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(fsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ssteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(gsteps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ststeps,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nrst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(sprst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(strst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!--------
#ifdef SOL_GPE
CALL MPI_BCAST(xi,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(csound,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(rho_avg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(disp_factor,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(oneybeta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nuN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(chempot0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#ifdef GPE_COUNTERFLOW
CALL MPI_BCAST(vnx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(vny,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef GPE_SELFTRUNC
CALL MPI_BCAST(ktrunc,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef GPE_FRCDISP
CALL MPI_BCAST(famp_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(kf1_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(kf1_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nhypvis1_gp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(nhypvis2_gp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(vishypo_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(vishype_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(friclowk_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef GPE_ATTRACT
CALL MPI_BCAST(gstrength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef GPPART
CALL MPI_BCAST(gp_N_obj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(gp_mass_part,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ppotd,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ppotV0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(ppotthick,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#ifdef GPOTOC
CALL MPI_BCAST(k_otoc_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BCAST(eps_otoc_gp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#endif
!--------

!! Set the physical and spectral space grid.
! That is define/initialize all the needed variables.
! E.g., dx,dy,dz, etc.
CALL set_grid

IF (myrank==0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! Grid parameters'
PRINT*, 'Box lenght:', 'Lx',lengthx,'Ly',lengthy
PRINT*, 'Spatial grid size'
PRINT*, 'Lx/Nx:',dx
PRINT*, 'Ly/Ny:',dy
PRINT*, '!!-----------------------------------------------!!'
PRINT*, 'Dealiasing: 2/3-rule (Spherical)'
PRINT*, 'Dealiasing wave number:', kalias
PRINT*, '!!-----------------------------------------------!!'
!! Data aquisition time
PRINT*,'Data acquisiton time (global, spectra, fields):', gtime,stime,ftime
ENDIF


!---------GPE-----------------------------------------

#ifdef SOL_GPE
!! Parameters characterising the GPE state.
Nparticle = lengthx*lengthy
xi = xi*dx
alpha = csound*xi*SQRT(rho_avg/two)
#ifndef GPE_ATTRACT
gstrength = csound**2/(two*alpha)
#endif
#ifndef GPE_SGLE
chempot0 = gstrength
omega_chem = chempot0
!IF (rho_avg .EQ. one) omega_chem = gstrength
#else
omega_chem = chempot0
#endif
!! 
noisesigma = SQRT((one/(two*alpha))*oneybeta)
!! Later a factor of sqrt(1/dvol) is multiplied to it in the SGLE.
!! Is it really a noisesigma then?

IF(myrank==0) PRINT*,gstrength,omega_chem,rho_avg,csound
IF(myrank==0) PRINT*,alpha,xi
IF(myrank==0) PRINT*,Nparticle

#ifdef GPE_FRCDISP
kf1_gp=kf1_gp*kdx
kf2_gp=kf2_gp*kdx
#endif

#ifdef GPPART
! Parameters of the potential modified.
ppotV0 = ppotV0*gstrength!*dexp(one/4.0_dbl)
ppotd = ppotd*xi
ppotthick = ppotthick*xi

CALL create_particles
#endif

#endif


!!-----------------------------------------------------

!! Data aquisition time
IF (myrank==0) PRINT*,'Data acquisiton time (global, spectra, fields):', gtime,stime,ftime


IF (myrank .EQ. 0) PRINT*, 'Starting FFTW based set up.'

!! Local array size estimation for mpi. (Contained in mod_fft.f90)
CALL mpi_localsize_estimate

IF (myrank .EQ. 0) PRINT*, 'Back in main after FFTW based size estimation.'


!!-----------------------------------------------------
!! Memory allocation and pointer association.
!------ State variables vel, omega, etc.
CALL create_state_variables

!! Create and set the global arrays, kx,ky, etc.
!--Create the arrays
CALL create_global_array

IF (myrank .EQ. 0) PRINT*, 'State variables allocated.'

!! Set up the Fourier transform operations (plans).
IF (myrank .EQ. 0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*,'Begin setting up FFTW plans.'
PRINT*, '!!-----------------------------------------------!!'
ENDIF
CALL initialize_fft

IF (myrank .EQ. 0) PRINT*,'FFT set-up complete. Back in main.'

!--Set the arrays
CALL global_array

!!-------------------------------------------------------
!! Set up the initial configuration, data, etc.         !
!!-------------------------------------------------------

IF (nrst .NE. 0) THEN
!!----------------------------------------
!! Case: Restart from an existing state  !
!!
!! nrst:: indicates the restart file no. !
!!----------------------------------------

IF (myrank == 0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*, 'Restart file number (nrst):', nrst
PRINT*, '!!-----------------------------------------------!!'
ENDIF

!!----------------------------------------
!! For Gross-Pitaevskii simulations.     !
!!----------------------------------------
#ifdef SOL_GPE
CALL restart_gpe(nrst)

!!----------------------------------------
!! GPE with particles                    !
!!----------------------------------------
#ifdef GPPART
CALL particle_potential
CALL initialize_particles
#endif

#endif

ELSE
!#else

IF (myrank == 0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! Starting at t=0'
PRINT*, 'Restart file number (nrst):', nrst
PRINT*, '!!-----------------------------------------------!!'
ENDIF

!!---------------------------------------
!! Initial wave functions for the GPE   !
!!---------------------------------------
#ifdef SOL_GPE
!CALL ic_gpe_TG
!CALL ic_gpe_kgaussrand
!CALL ic_gpe
!CALL ic_gpe_vortex
!CALL ic_combine2wf

!!----------------------------------------
!! If want to restart from ARGLE or SGLE !
!! prepared data.                        !
!!----------------------------------------
CALL restart_from_argle_data(nrst)
#ifdef GPOTOC
CALL ic_otoc_gpe
#endif

!!----------------------------------------
!! GPE with particles.                   !
!!----------------------------------------
#ifdef GPPART
!! (i) Set the particle potential
CALL particle_potential
!! (ii) Position, velocity and force.
CALL initialize_particles
#endif

#endif
ENDIF

!! Solver: contains time iterations, io's etc.
CALL solve

!! Destroy the plans created for the Fourier transfrom
!! operations.
CALL clear_fft

!! Release the allocated memory.
CALL free_state_variables
CALL dealloc_global_array

#ifdef GPPART
CALL destroy_particles
#endif

!=====================================================================
        
! Calculation of time taken
CALL mpi_barrier(MPI_COMM_WORLD,ierr)
time_fin = mpi_wtime()

IF (myrank == 0) THEN
OPEN(UNIT=1,FILE='run_info.dat',STATUS='unknown',IOSTAT=ios)
WRITE(1,*) 'time taken by',myrank,'processor',&
time_fin-time_ini
CLOSE(1)
ENDIF

call MPI_FINALIZE(ierr)

END PROGRAM
