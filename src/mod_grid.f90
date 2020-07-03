MODULE mod_grid

USE mod_precision
USE mod_constants
USE FFTW3 !! Needed for fftw data types.

IMPLICIT NONE

!! Glossary

CHARACTER(100)::fnum,prc
CHARACTER(100)::fnrst

INTEGER:: ios

INTEGER, PARAMETER:: ndim=2

INTEGER:: Nx
INTEGER:: Ny
!INTEGER(C_INTPTR_T):: Nx
!INTEGER(C_INTPTR_T):: Ny
!INTEGER(C_INT):: Nx
!INTEGER(C_INT):: Ny

!INTEGER(C_INTPTR_T):: Nxh,Nyh,Nzh
!INTEGER(C_INTPTR_T):: Nxhp,Nyhp,Nzhp
!INTEGER(C_INTPTR_T):: Nxpp,Nypp,Nzpp

!! No. of files into which data is split at run time originally.
INTEGER:: nsplit_orignal,nsplit_current

INTEGER:: Nxh,Nyh
INTEGER:: Nxhp,Nyhp
INTEGER:: Nxpp,Nypp

INTEGER:: ksqr_max,nshell

REAL(KIND=GP):: lengthx,lengthy
REAL(KIND=GP):: dx,dy
REAL(KIND=GP):: kdx,kdy
REAL(KIND=GP):: dvol
REAL(KIND=GP):: kdvol
REAL(KIND=GP):: kalias,kasqr

INTEGER:: maxiter
INTEGER:: fsteps,ssteps,gsteps,ststeps
INTEGER:: nrst,sprst,strst
INTEGER:: tstart
REAL(KIND=GP):: ftime,stime,gtime
REAL(KIND=GP):: timestart
REAL(KIND=GP):: dt
!! Fixed time step.
REAL(KIND=GP):: dt_fix
!! Adoptive time step computed for Navier-Stokes fluid.
REAL(KIND=GP):: dt_adpt_ns 
REAL(KIND=GP):: dt_adpt_2f
REAL(KIND=GP):: dt_adpt_scl1
!! CFL tolerance.
REAL(KIND=GP):: cflmax

CONTAINS

SUBROUTINE set_grid

IMPLICIT NONE

Nxh = Nx/2
Nyh = Ny/2

Nxhp = Nx/2+1
Nyhp = Ny/2+1
Nxpp = Nx+2
Nypp = Ny+2

ksqr_max = (Nxh**2+Nyh**2)

nshell = int(1.414*Nx/2.0d0) + 1

lengthx = two*pi
lengthy = two*pi

dx = lengthx/Nx
dy = lengthy/Ny

kdx = two*pi/lengthx
kdy = two*pi/lengthy

dvol = dx*dy
kdvol = kdx*kdy

! Assume Nx=Ny
kalias = REAL(Nx,KIND=GP)/3.0_GP
!kalias=12

!! Set the data aquisition time 
gtime = gsteps*dt_fix
stime = ssteps*dt_fix
ftime = fsteps*dt_fix


END SUBROUTINE set_grid

END MODULE mod_grid
