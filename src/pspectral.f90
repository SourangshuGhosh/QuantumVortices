SUBROUTINE global_array

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD
USE mod_mpi

IMPLICIT NONE

INTEGER:: i1,i2,i3,i4
INTEGER:: i2_loc,i3_loc
INTEGER:: ishell
INTEGER:: k1,k2,k3,ksqr,mshl
REAL(KIND=GP):: kx,ky,kz
REAL(KIND=GP)::rk,rk2


#ifdef SOL_GPE
DO i1 = 1,Nx
! FT complex to complex
kwn1_c2c(i1) = (i1-1) - Nx*(i1/Nxhp)
!if (myrank ==0) print*,'k1',kwn1_c2c(i1)
ENDDO

DO i2 = 1,Ny
kwn2_c2c(i2) = (i2-1) - Ny*(i2/Nyhp)
!if (myrank ==0) print*,'k2',kwn2_c2c(i2)
ENDDO

#endif

IF (myrank==0) PRINT*,'Global array subroutine is okay'

END SUBROUTINE global_array



subroutine init_random_seed(iopoint,myrank)

use iso_fortran_env, only: int64
implicit none
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
integer(int64) :: t

CHARACTER(LEN=100):: prc,fnum
INTEGER:: i1,ios
INTEGER, INTENT(IN):: myrank,iopoint
REAL*8 :: rrn(1:20)
          
call random_seed(size = n)
            
allocate(seed(n))
            
! First try if the OS provides a random number generator
!            open(newunit=un, file="/dev/urandom", access="stream", &
!                 form="unformatted", action="read", status="old", iostat=istat)
!            if (istat == 0) then
!               read(un) seed
!               close(un)
!            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
!               call system_clock(t)
!               print*,myrank,t
!               if (t == 0) then
!                  call date_and_time(values=dt)
!                  print*,myrank, dt
!                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
!                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
!                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
!                       + dt(5) * 60 * 60 * 1000 &
!                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
!                       + dt(8)
!               end if
!               !pid = getpid()
!               !t = ieor(t, int(pid, kind(t)))
!               do i = 1, n
!                  seed(i) = lcg(t)
!               end do
!!            end if

!! I would prefer to use fixed date.
!CALL DATE_AND_TIME(VALUES=dt)
dt(1) = 2020
dt(2) = 4
dt(3) = 18
dt(4) = 330
dt(5) = 20
dt(6) = 30
dt(7) = 56
dt(8) = 468
t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
     + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
     + dt(3) * 24_int64 * 60 * 60 * 1000 &
     + dt(5) * 60 * 60 * 1000 &
     + dt(6) * 60 * 1000 + dt(7) * 1000 &
     + dt(8)+ myrank
PRINT*,myrank,t
DO i = 1, n
seed(i) = lcg(t)
ENDDO

WRITE(fnum,'(i8)')iopoint 
WRITE(prc,'(i8)') myrank+1

OPEN(UNIT=1,FILE='data/ranf_seeds/seed_ini_'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED',STATUS='NEW',IOSTAT=ios)
DO i1=1,n
WRITE(1) seed(i1)
ENDDO
CLOSE(1)

!OPEN(UNIT=1,FILE='my_seed_p'//TRIM(ADJUSTL(prc))//'.dat',STATUS='NEW',IOSTAT=ios)
!DO i1=1,n
!WRITE(1,*) seed(i1)
!ENDDO
!CLOSE(1)               


contains
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
function lcg(s)
integer :: lcg
integer(int64) :: s
if (s == 0) then
s = 104729
else
s = mod(s, 4294967296_int64)
end if
s = mod(s * 279470273_int64, 4294967291_int64)
lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg
end subroutine init_random_seed

