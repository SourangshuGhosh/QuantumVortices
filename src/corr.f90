PROGRAM corr

USE FFTW3

IMPLICIT NONE

INCLUDE "mpif.h"

! Some constants
integer, parameter:: GP = KIND(0.0D0)
real(kind=GP), parameter:: zero  = 0.0_GP
real(kind=GP), parameter:: one   = 1.0_GP
real(kind=GP), parameter:: mone  = -1.0_GP
real(kind=GP), parameter:: pi    = 4.0_GP*atan(1.0_GP)
real(kind=GP), parameter:: half  = 1.0_GP/2.0_GP
real(kind=GP), parameter:: oney4 = 1.0_GP/4.0_GP
real(kind=GP), parameter:: two   = 2.0_GP
real(kind=GP), parameter:: oney3 = 1.0_GP/3.0_GP 
real(kind=GP), parameter:: oney6 = 1.0_GP/6.0_GP
complex(kind=GP), parameter:: czero = CMPLX(0.0,0.0,KIND=GP)!complex(0.0_GP,0.0_GP)
complex(kind=GP), parameter:: zi    = CMPLX(0.0,1.0,KIND=GP)!complex(0.0_GP,1.0_GP)

!! MPI
!-------------------------------------------------------
INTEGER:: istatus(MPI_STATUS_SIZE)
INTEGER:: ierr,nprocs,myrank
INTEGER:: time_ini, time_fin
!--------------------------------------------------------

INTEGER(KIND=8) :: rec_len
INTEGER:: Nx,Ny,Nz,Nt,Ntt
INTEGER:: Nxh,Nyh,Nzh
INTEGER:: Nxhp,Nyhp,Nzhp
INTEGER:: Nxpp,Nypp,Nzpp
INTEGER:: mshl,nshell
INTEGER:: ios
INTEGER:: i1,i2,i3,ic1,ic2,ic3,cstep
INTEGER:: i4
INTEGER:: ix,iy,iz,ir
INTEGER:: k1,k2,k3
INTEGER:: ksqr
INTEGER:: factor1,factor2
INTEGER:: ifile
INTEGER:: cloop,cntfile
INTEGER:: cxmin,cxmax,cymin,cymax,czmin,czmax
INTEGER:: ip, order

REAL(KIND=GP):: lengthx,lengthy,lengthz
REAL(KIND=GP):: x,y,z
REAL(KIND=GP):: dx,dy,dz
REAL(KIND=GP):: du1,du2,du3
REAL(KIND=GP):: su
REAL(KIND=GP):: rs
REAL(KIND=GP):: k0,a,b,c

REAL(KIND=GP):: ksqrt

INTEGER(C_INTPTR_T):: RNx
INTEGER(C_INTPTR_T):: RNy

INTEGER(C_INTPTR_T):: alloc_local1, alloc_local3,locsize
INTEGER(C_INTPTR_T):: Rlocal_Ny,Rlocal_i2_offset
INTEGER(C_INTPTR_T):: Rlocal_Nx,Rlocal_i1_offset
INTEGER(C_INTPTR_T):: nhowmany
INTEGER(C_INTPTR_T), PARAMETER:: dim2d=2
INTEGER(C_INTPTR_T):: dimnd, nscl,nv
! For complex to complex transforms
INTEGER(C_INTPTR_T):: alloc_local_cmpx
INTEGER:: local_Ny_cmpx, local_i2_offset_cmpx
INTEGER:: local_Nx_cmpx, local_i1_offset_cmpx

!
TYPE::plan
CHARACTER (100):: pname
INTEGER:: pflag
TYPE(C_PTR):: pfor
TYPE(C_PTR):: pinv
INTEGER(C_INT)::rank
INTEGER(C_INTPTR_T), DIMENSION(1:2)::nrank
END TYPE
! Local size computed using FFTW subroutine.
TYPE:: lclsize
INTEGER(C_INT):: rank
INTEGER(C_INTPTR_T), DIMENSION(1:2):: nrank
END TYPE

INTEGER, PARAMETER:: forward=-1, inverse=1
REAL(KIND=GP):: fftwscale

TYPE(lclsize):: lclsize2

TYPE(plan):: plan_psi

TYPE(C_PTR):: data_psi,data_kpsi
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: psi(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: kpsi(:,:)

!========================================================
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

Nx = 256
Ny = 256
Nz = 256

Nxh = Nx/2
Nyh = Ny/2
Nzh = Nz/2

Nxhp = Nx/2+1
Nyhp = Ny/2+1
Nzhp = Nz/2+1
Nxpp = Nx+2
Nypp = Ny+2
Nzpp = Nz+2

! Box length
lengthx = two*pi
lengthy = two*pi
lengthz = two*pi

dx = lengthx/Nx
dy = lengthy/Ny
dz = lengthz/Nz

!=============================================================

!! Obtain the size of the local arrays on each process, using FFTW
!! slab decomposition.
!!--------------------------------------------
RNx=Nx
RNy=Ny

lclsize2%rank = 2
lclsize2%nrank(1:2) = (/RNy,RNx/)

alloc_local_cmpx = fftw_mpi_local_size_transposed(lclsize2%rank,lclsize2%nrank,&
MPI_COMM_WORLD,Rlocal_Ny,Rlocal_i2_offset,Rlocal_Nx,Rlocal_i1_offset)
local_Ny_cmpx=Rlocal_Ny
local_Nx_cmpx=Rlocal_Nx
local_i2_offset_cmpx=Rlocal_i2_offset
local_i1_offset_cmpx=Rlocal_i1_offset
!------------------------------------------------

data_psi = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_psi,psi,[Nx,local_Ny_cmpx])

data_kpsi = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_kpsi,kpsi,[Ny,local_Nx_cmpx])
!------------------------------------------------

!! FFTW plans
!------------------------------------------------
plan_psi%pname = 'psi'
plan_psi%pflag = 1
plan_psi%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_psi%nrank(1:2)=(/RNy,RNx/)

IF (myrank .EQ. 0) PRINT*, 'Plan psi'

!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_psi%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       psi,kpsi,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_psi%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       kpsi,psi,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)


fftwscale = one/(REAL(Nx,KIND=GP)*REAL(Ny,KIND=GP))
!---------------------------------------------------


!CALL fftw_mpi_execute_dft(plan_psi%pfor,psi,kpsi)
!ELSEIF (plan_name .EQ. plan_psi%pname .AND. direction .EQ. inverse) THEN
!CALL fftw_mpi_execute_dft(plan_psi%pinv,kpsi,psi)


!! Destroy the plans
!-------------------------------------
CALL fftw_destroy_plan(plan_psi%pfor)
CALL fftw_destroy_plan(plan_psi%pinv)

!-------------------------------------
CALL fftw_free(data_psi)
CALL fftw_free(data_kpsi)


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

CALL MPI_FINALIZE(ierr)

END PROGRAM 
