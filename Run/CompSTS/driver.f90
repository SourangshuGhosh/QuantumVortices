PROGRAM driver

USE FFTW3

IMPLICIT NONE

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

COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: sts1d1
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: arr1

REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: spec_sts2d

REAL(KIND=GP):: enorm,stsnorm

INTEGER:: r_iopoint
REAL(KIND=GP):: r_time,r_dt

CHARACTER(100)::fnn,prcjj
INTEGER:: filen,nregrpfile,ii,jj,kk,nplnperproc_orignal
INTEGER:: nsplit_original,nsplit_current
INTEGER:: i1_loc,i2_loc,myrank,local_Nx_cmpx


!! Hann function.
REAL(KIND=GP):: whann, tpp

!! For FFTW.
INTEGER(C_INT):: RNt,RNthm
TYPE(C_PTR):: pfor_sts1d
COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:):: ft1d1
REAL(KIND=GP):: fftwscale

INTEGER(C_INT):: rank
INTEGER(C_INT), DIMENSION(1:1):: nrank
INTEGER(C_INT):: howmany
INTEGER(C_INT), DIMENSION(1:1):: inembed,onembed
INTEGER(C_INT):: istride, ostride
INTEGER(C_INT):: idist, odist
TYPE(C_PTR):: pfor_sts2dhm
COMPLEX(C_DOUBLE_COMPLEX), ALLOCATABLE, DIMENSION(:,:,:):: sts2d1

! -- idist, odist
! The distance in memory between the first element of the first array 
! and the first element of the second array */

Nx = 128
Ny = 128
Nt = 2048
Ntt = 2048
RNt = 2048
RNthm = Ntt


Nxh = Nx/2
Nyh = Ny/2

Nxhp = Nx/2+1
Nyhp = Ny/2+1
Nxpp = Nx+2
Nypp = Ny+2

nshell = int(1.4142*Nx/2.0d0) + 1

lengthx = two*pi
lengthy = two*pi

dx = lengthx/Nx
dy = lengthy/Ny


nsplit_original = 4
nsplit_current = 1 
myrank = 0
filen = 50

local_Nx_cmpx = Nx/nsplit_original

ALLOCATE(sts1d1(1:Nx,1:Nt))
ALLOCATE(arr1(1:Nx,1:Ntt))
ALLOCATE(ft1d1(1:Ntt))

ALLOCATE(sts2d1(1:Nx,1:Ny,1:Nt))
ALLOCATE(spec_sts2d(0:nshell,1:Ntt))

enorm = (one/REAL(Nx*Ny,KIND=GP))**2
fftwscale = one/real(Ntt,kind=GP)
stsnorm = (one/REAL(Ntt,KIND=GP))**2*(one/REAL(Nx*Ny,KIND=GP))**2
print*,Ntt,stsnorm


!! FFTW 1D

pfor_sts1d = fftw_plan_dft_1d(RNt,ft1d1,ft1d1,FFTW_FORWARD,FFTW_MEASURE)

!! FFTW 1D Howmany
rank = 1
nrank(1:1) = (/Ntt/)
howmany = Nx*Ny
idist = 1
odist = 1
istride = Nx*Ny
ostride = Nx*Ny
inembed = nrank
onembed = nrank
pfor_sts2dhm = fftw_plan_many_dft(rank,nrank,howmany,sts2d1,inembed,istride,idist,&
                                sts2d1,onembed,ostride,odist,FFTW_FORWARD,FFTW_MEASURE)


!! Read the kpsi kx=0 

WRITE(fnn,'(i8)') filen
print*,fnn

OPEN(UNIT=1,FILE='../data/sts1d/kpsi_x0_nrst'//TRIM(ADJUSTL(fnn))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

DO i4 = 1,Nt
READ(1) (sts1d1(i2,i4),i2=1,Ny)
ENDDO

CLOSE(1)

!!--(kx=0)
DO i2 = 1,Ny
DO i4 = 1,Ntt
tpp = 2*pi*(i4-1)/(Ntt-1)
whann = half*(one - COS(tpp) )
ft1d1(i4) = whann*sts1d1(i2,i4)
!print*, abs(ft1d1(i4))
ENDDO

CALL FFTW_EXECUTE_DFT(pfor_sts1d,ft1d1,ft1d1)

DO i4 = 1,Ntt
arr1(i2,i4) = stsnorm*ABS(ft1d1(i4))**2
!write(i1+10,*), abs(ft1d1(i4))
ENDDO

ENDDO


!INQUIRE (IOLENGTH=rec_len) arr1
!!--(kx,0,0)
!OPEN(UNIT=1,FILE='sts1d_kx0_nrst'//TRIM(ADJUSTL(fnn))//'.dat',&
!FORM='UNFORMATTED', status='NEW', IOSTAT=ios, ACCESS='direct',RECL=rec_len)
!!DO i4 = 1,Ntt
!!WRITE(1,rec=1) (arr1(i1,i4),i1=1,Nx)
!WRITE(1,rec=1) arr1
!!ENDDO
!CLOSE(1)

OPEN(UNIT=1,FILE='sts1d_kx0_nrst'//TRIM(ADJUSTL(fnn))//'.dat',&
 status='NEW', IOSTAT=ios)
DO i4 = 1,Ntt
WRITE(1,*) (arr1(i2,i4),i2=1,Ny)
ENDDO
CLOSE(1)


!----------------------------------

!! Planes

WRITE(fnn,'(i8)') filen
print*,fnn

!--- Read the data
DO i3 = 1,nsplit_original

WRITE(prcjj,'(i8)') i3
print*,prcjj

OPEN(UNIT=1,FILE='../data/sts2d/kpsi_nrst'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

DO i4 = 1,Nt
READ(1) ((sts2d1(i2,i1,i4),i2=1,Ny),i1=(i3-1)*local_Nx_cmpx+1,(i3)*local_Nx_cmpx)
ENDDO

CLOSE(1)
ENDDO
!--------

DO i4 = 1,Ntt
DO i1 = 1,Nx
k1 = (i1-1) - Nx*(i1/Nxhp)
DO i2 = 1,Ny
k2 = (i2-1) - Ny*(i2/Nyhp)

tpp = 2*pi*(i4-1)/(Ntt-1)
whann = half*(one - COS(tpp) )
sts2d1(i2,i1,i4) = whann*sts2d1(i2,i1,i4)
ENDDO
ENDDO
ENDDO

CALL FFTW_EXECUTE_DFT(pfor_sts2dhm,sts2d1,sts2d1)

spec_sts2d = zero

DO i4 = 1,Ntt
DO i1 = 1,Nx
k1 = (i1-1) - Nx*(i1/Nxhp)
DO i2 = 1,Ny
k2 = (i2-1) - Ny*(i2/Nyhp)

ksqr = k1**2+k2**2
ksqrt = SQRT(REAL(ksqr,KIND=GP))
mshl = NINT(ksqrt)

!arr3(i1,i4) =  stsnorm*ABS(sts2d1(i1,i3,i4))**2

spec_sts2d(mshl,i4) = spec_sts2d(mshl,i4) + stsnorm*ABS(sts2d1(i2,i1,i4))**2

ENDDO
ENDDO
ENDDO

OPEN(UNIT=1,FILE='sts2d_nrst'//TRIM(ADJUSTL(fnn))//'.dat',&
status='NEW', IOSTAT=ios)
DO i4 = 1,Ntt
WRITE(1,*) (spec_sts2d(mshl,i4),mshl=0,INT(Nx/3))
ENDDO
CLOSE(1)

CALL FFTW_DESTROY_PLAN(pfor_sts1d)
CALL FFTW_DESTROY_PLAN(pfor_sts2dhm)

DEALLOCATE(sts1d1)
DEALLOCATE(arr1)
DEALLOCATE(ft1d1)
DEALLOCATE(sts2d1)
DEALLOCATE(spec_sts2d)

!======================================================================

END PROGRAM
