PROGRAM rho_vzu

IMPLICIT NONE

INTEGER(KIND=8) :: rec_len
INTEGER:: Nx,Ny,Nz
INTEGER:: Nxh,Nyh,Nzh
INTEGER:: Nxhp,Nyhp,Nzhp
INTEGER:: Nxpp,Nypp,Nzpp
INTEGER:: nshell
INTEGER:: ios
INTEGER:: i1,i2,i3,ic1,ic2,ic3,cstep
INTEGER:: ix,iy,iz,ir
INTEGER:: factor1,factor2
INTEGER:: ifile
INTEGER:: cloop,cntfile
INTEGER:: cxmin,cxmax,cymin,cymax,czmin,czmax
INTEGER:: ip, order
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

REAL(KIND=GP):: lengthx,lengthy,lengthz
REAL(KIND=GP):: x,y,z
REAL(KIND=GP):: dx,dy,dz
REAL(KIND=GP):: du1,du2,du3
REAL(KIND=GP):: su
REAL(KIND=GP):: rs
REAL(KIND=GP):: k0,a,b,c
COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: psi
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: rho

INTEGER:: r_iopoint
REAL(KIND=GP):: r_time,r_dt

CHARACTER(100)::fnn,prcjj
INTEGER:: filen,nregrpfile,ii,jj,kk,nplnperproc_orignal
INTEGER:: nsplit_orignal,nsplit_current
INTEGER:: i2_loc,i3_loc,myrank

Nx = 256
Ny = 256

Nxh = Nx/2
Nyh = Ny/2

Nxhp = Nx/2+1
Nyhp = Ny/2+1
Nxpp = Nx+2
Nypp = Ny+2

lengthx = two*pi
lengthy = two*pi

dx = lengthx/Nx
dy = lengthy/Ny


nsplit_orignal = 4
nsplit_current = 1 
myrank = 0
filen = 8

ALLOCATE(psi(1:Nx,1:Ny))
ALLOCATE(rho(1:Nx,1:Ny))

OPEN(UNIT=12,FILE='rho/rho'//TRIM(ADJUSTL(fnn))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios,ACCESS='direct')
READ(12) rho
CLOSE(12)


!do i1=1,Nx
!write(11,*) rho(i1,32,32)
!enddo

!INQUIRE (IOLENGTH=rec_len) rho


DEALLOCATE(psi)
DEALLOCATE(rho)

END PROGRAM
