#ifdef SOL_GPE

SUBROUTINE ic_gpe_kgaussrand

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_fft
USE mod_2dHD
USE mod_mpi
USE mod_grid
USE mtmod

IMPLICIT NONE
INTEGER:: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: mshl,nn
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius,kradius2
REAL(KIND=GP):: k0
REAL(KIND=GP):: norm,tempnorm,alphak,ran,stndev,prefac

stndev = half*kdx
prefac = one/(dsqrt(dsqrt(pi)*stndev))
nn = 2
k0 = nn*kdx      !! This is an adjustable parameter.

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

ksqr = k1**2+k2**2
kradius2 = kx**2 + ky**2
kradius = sqrt(kradius2)

CALL RANDOM_NUMBER(ran)
!ran = grnd()
alphak = ran*two*pi
kpsi(i2,i1_loc) = prefac*exp(-((kradius-k0)**2)/(two*stndev**2))*exp(zi*alphak)

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

CALL fft(plan_psi%pname,inverse)

! Normalization of the generated wavefunction.
norm = zero
DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
norm = norm + (abs(psi(i1,i2_loc)))**2
ENDDO
ENDDO

CALL MPI_REDUCE(norm,tempnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) norm = tempnorm*dvol
CALL MPI_BCAST(norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = psi(i1,i2_loc)*sqrt(lengthx*lengthy)/sqrt(norm)
ENDDO
ENDDO

CALL fft(plan_psi%pname,forward)

END SUBROUTINE ic_gpe_kgaussrand

SUBROUTINE ic_gpe

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_fft
USE mod_2dHD
USE mod_mpi
USE mod_grid

IMPLICIT NONE
INTEGER:: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: mshl,nn
REAL(KIND=GP):: x,y
REAL(KIND=GP):: norm,tempnorm
REAL(KIND=GP):: ax,ay,tterm

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
y = (i2-1)*dy
DO i1 = 1, Nx
x = (i1-1)*dx

!az = COS(z*sqrt(two))
!ay = COS(y*sqrt(two))
!tterm = az + zi*ay
!psi(i1,i2,i3_loc) = tterm*TANH(ABS(tterm)/(sqrt(two)*xi))/ABS(tterm)

psi(i1,i2_loc) = one!+1.0d-2*COS(x)**2!exp(zi*(i1-1)*dx)

!if (x .GE. pi+15*dx .and. x .LE. pi+30*dx .and. y .GE. pi+5*dx .and. y .LE. pi+15*dx) psi(i1,i2_loc) = czero

ENDDO
ENDDO

CALL fft(plan_psi%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

CALL fft(plan_psi%pname,inverse)

! Normalization of the generated wavefunction.
norm = zero
DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
norm = norm + (abs(psi(i1,i2_loc)))**2
ENDDO
ENDDO

CALL MPI_REDUCE(norm,tempnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) norm = tempnorm*dvol
CALL MPI_BCAST(norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = psi(i1,i2_loc)*sqrt(lengthx*lengthy)/sqrt(norm)
ENDDO
ENDDO

CALL fft(plan_psi%pname,forward)

END SUBROUTINE ic_gpe

SUBROUTINE ic_gpe_TG

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_fft
USE mod_2dHD
USE mod_mpi
USE mod_grid

IMPLICIT NONE
INTEGER:: time
INTEGER:: i1,i2,i3
INTEGER:: i1_loc,i2_loc,i3_loc
INTEGER:: k1,k2,k3,ksqr
INTEGER:: mshl,nn
REAL(KIND=GP):: x,y
REAL(KIND=GP):: norm,tempnorm
REAL(KIND=GP):: isum, tsum, psign
REAL(KIND=GP):: gammad
REAL(KIND=GP):: cleb1, cleb2
REAL(KIND=GP):: ampcleb, oneysqrt2
COMPLEX(KIND=GP):: term1,term2,term3,term4

gammad = int(8.0d0/(4.0*(4.0d0*pi*alpha)))
IF (myrank .EQ. 0) PRINT*, 'gammad',gammad
oneysqrt2 = dsqrt(one/two)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
y = (i2-1)*dy
DO i1 = 1, Nx
x = (i1-1)*dx

        cleb1 = COS(x)*SQRT(two)
        cleb2 = COS(y)*SQRT(two)

        ampcleb = SQRT((cleb1-oneysqrt2)**2+cleb2**2)
        term1 = dcmplx(cleb1-oneysqrt2,cleb2)*&
                        TANH(ampcleb/(SQRT(two)*xi))/ampcleb

        ampcleb = SQRT(cleb1**2+(cleb2-oneysqrt2)**2)
        term2 = dcmplx(cleb1,cleb2-oneysqrt2)*&
                        TANH(ampcleb/(SQRT(two)*xi))/ampcleb

        ampcleb = SQRT((cleb1+oneysqrt2)**2+cleb2**2)
        term3 = dcmplx(cleb1+oneysqrt2,cleb2)*&
                        TANH(ampcleb/(SQRT(two)*xi))/ampcleb

        ampcleb = SQRT(cleb1**2+(cleb2+oneysqrt2)**2)
        term4 = dcmplx(cleb1,cleb2+oneysqrt2)*&
                        TANH(ampcleb/(SQRT(two)*xi))/ampcleb

        psi(i1,i2_loc) = term1*term2*term3*term4

ENDDO
ENDDO

CALL fft(plan_psi%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

CALL fft(plan_psi%pname,inverse)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx

psi(i1,i2_loc) = psi(i1,i2_loc)**gammad

ENDDO
ENDDO

CALL fft(plan_psi%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

CALL fft(plan_psi%pname,inverse)

! Normalization of the generated wavefunction.
norm = zero
DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
norm = norm + (abs(psi(i1,i2_loc)))**2
ENDDO
ENDDO

CALL MPI_REDUCE(norm,tempnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) norm = tempnorm*dvol
CALL MPI_BCAST(norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = psi(i1,i2_loc)*sqrt(lengthx*lengthy)/sqrt(norm)
ENDDO
ENDDO

CALL fft(plan_psi%pname,forward)

END SUBROUTINE ic_gpe_TG

#ifdef GPOTOC

SUBROUTINE ic_OTOC_gpe

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_fft
USE mod_2dHD
USE mod_mpi
USE mod_grid
USE mtmod

IMPLICIT NONE
INTEGER:: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: mshl,nn
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius,kradius2
REAL(KIND=GP):: k0
REAL(KIND=GP):: norm,tempnorm

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

ksqr = k1**2+k2**2
rk2  = kx**2+ky**2
rk = sqrt(rk2)
mshl = nint(rk)

kradius2 = kx**2 + ky**2

IF (mshl == k_otoc_gp) THEN
kpsi(i2,i1_loc) = kpsi(i2,i1_loc)*(one+eps_otoc_gp)
ENDIF

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

CALL fft(plan_psi%pname,inverse)

! Normalization of the generated wavefunction.
!norm = zero
!DO i2_loc = 1, local_Ny_cmpx
!i2 = i2_loc + local_i2_offset_cmpx
!DO i1 = 1, Nx
!norm = norm + (abs(psi(i1,i2_loc)))**2
!ENDDO
!ENDDO

!CALL MPI_REDUCE(norm,tempnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
!                0,MPI_COMM_WORLD,ierr)
!IF (myrank .EQ. 0) norm = tempnorm*dvol
!CALL MPI_BCAST(norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!DO i2_loc = 1, local_Ny_cmpx
!i2 = i2_loc + local_i2_offset_cmpx
!DO i1 = 1, Nx
!psi(i1,i2_loc) = psi(i1,i2_loc)*sqrt(lengthx*lengthy)/sqrt(norm)
!ENDDO
!ENDDO

!CALL fft(plan_psi%pname,forward)

END SUBROUTINE ic_OTOC_gpe

#endif

SUBROUTINE ic_combine2wf

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi

IMPLICIT NONE

INTEGER:: i1,i2,i3,i4
INTEGER:: i1_loc,i2_loc
INTEGER:: ishell
INTEGER:: k1,k2,ksqr,mshl
!INTEGER:: time,icpoint
!INTEGER, INTENT(IN):: nrestart
INTEGER:: icpoint
INTEGER:: jtime
REAL(KIND=GP):: time
REAL(KIND=GP):: dt_argle
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky
REAL(KIND=GP):: rk,rk2
REAL(KIND=GP):: atemp,amin,amax,tmax,tmin
REAL(KIND=GP):: u1,u2,u3,u4,u5,u6,ss1
COMPLEX(KIND=GP):: cwf
REAL(KIND=GP):: norm,tempnorm
!! Regrouping of data files.
CHARACTER(100)::fnn,prcjj
INTEGER:: nregrpfile,ii,jj,kk,nplnperproc_orignal

CHARACTER(LEN=100):: filename1,filename2

!------------------------------------------------------------------------------------
#ifdef IOEVERYPROC
! Specify partial names: "pX.dat" part would be completed in read_field subroutine.
filename1 = 'ARGLE/wfc1'
CALL read_field(filename1)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = cgp_phys_tmp1(i1,i2_loc)
ENDDO
ENDDO

filename2 = 'ARGLE/wfc2'
CALL read_field(filename2)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = psi(i1,i2_loc)*cgp_phys_tmp1(i1,i2_loc)
ENDDO
ENDDO
#endif

#ifdef IOHDF5
! Specify full name.
filename1 = 'ARGLE/wfc1.h5'
CALL read_field(filename1)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = cgp_phys_tmp1(i1,i2_loc)
ENDDO
ENDDO

filename2 = 'ARGLE/wfc2.h5'
CALL read_field(filename2)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = psi(i1,i2_loc)*cgp_phys_tmp1(i1,i2_loc)
ENDDO
ENDDO
#endif

CALL fft(plan_psi%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

CALL fft(plan_psi%pname,inverse)

! Normalization of the generated wavefunction.
!norm = zero
!DO i2_loc = 1, local_Ny_cmpx
!i2 = i2_loc + local_i2_offset_cmpx
!DO i1 = 1, Nx
!norm = norm + (abs(psi(i1,i2_loc)))**2
!ENDDO
!ENDDO

!CALL MPI_REDUCE(norm,tempnorm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
!                0,MPI_COMM_WORLD,ierr)
!IF (myrank .EQ. 0) norm = tempnorm*dvol
!CALL MPI_BCAST(norm, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

!DO i2_loc = 1, local_Ny_cmpx
!i2 = i2_loc + local_i2_offset_cmpx
!DO i1 = 1, Nx
!psi(i1,i2_loc) = psi(i1,i2_loc)*sqrt(lengthx*lengthy)/sqrt(norm)
!ENDDO
!ENDDO

!CALL fft(plan_psi%pname,forward)

IF (myrank .EQ. 0) THEN
PRINT*,'Initial WF formed by combining: wfc1 and wfc2'
ENDIF

!! Set the time and file counter to zero.
timestart = zero
nrst = 0

END SUBROUTINE ic_combine2wf

SUBROUTINE read_field(filename)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi

IMPLICIT NONE

INTEGER:: i1,i2,i3,i4
INTEGER:: i1_loc,i2_loc
INTEGER:: ishell
INTEGER:: k1,k2,ksqr,mshl
!INTEGER:: time,icpoint
!INTEGER, INTENT(IN):: nrestart
INTEGER:: icpoint
INTEGER:: jtime
REAL(KIND=GP):: time
REAL(KIND=GP):: dt_argle
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky
REAL(KIND=GP):: rk,rk2
REAL(KIND=GP):: atemp,amin,amax,tmax,tmin
REAL(KIND=GP):: u1,u2,u3,u4,u5,u6,ss1
COMPLEX(KIND=GP):: cwf
!! Regrouping of data files.
CHARACTER(100)::fnn,prcjj
INTEGER:: nregrpfile,ii,jj,kk,nplnperproc_orignal

CHARACTER(LEN=100), INTENT(IN):: filename        

#ifdef IOEVERYPROC

IF (nsplit_orignal .NE. nsplit_current .AND. nsplit_orignal .GT. nsplit_current) THEN

nregrpfile=nsplit_orignal/nsplit_current
nplnperproc_orignal = Ny/nsplit_orignal

i2_loc=0

DO ii = 1,nregrpfile
jj = nregrpfile*myrank + ii

WRITE(prcjj,'(i8)') jj

OPEN(UNIT=11,FILE=TRIM(ADJUSTL(filename))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) icpoint,time,dt_argle

DO kk = 1,nplnperproc_orignal
i2_loc = i2_loc + 1

READ(11) (cgp_phys_tmp1(i1,i2_loc),i1=1,Nx)

ENDDO

CLOSE(11)

ENDDO

ELSEIF (nsplit_orignal .NE. nsplit_current .AND. nsplit_orignal .LT. nsplit_current) THEN
 
jj = ((Ny/nsplit_current)*myrank)/(Ny/nsplit_orignal)+1
WRITE(prcjj,'(i8)') jj

OPEN(UNIT=11,FILE=TRIM(ADJUSTL(filename))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) icpoint,time,dt_argle

i2_loc=0
DO ii = (jj-1)*(Ny/nsplit_orignal)+1, jj*(Ny/nsplit_orignal)

IF (ii .GT. (Ny/nsplit_current)*myrank .AND. ii .LE. (Ny/nsplit_current)*(myrank+1) ) THEN
i2_loc = i2_loc + 1

READ(11) (cgp_phys_tmp1(i1,i2_loc),i1=1,Nx)

ELSE

READ(11) (cwf,i1=1,Nx)

ENDIF
ENDDO

CLOSE(11)

ELSE

WRITE(prc,'(i8)') (myrank+1)

OPEN(UNIT=11,FILE=TRIM(ADJUSTL(filename))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) icpoint,time,dt_argle

DO i2_loc=1,local_Ny_cmpx
READ(11) (cgp_phys_tmp1(i1,i2_loc),i1=1,Nx)
ENDDO
CLOSE(11)

ENDIF

#endif

#ifdef IOHDF5
!filename = 'restart/wf/wf_'//TRIM(ADJUSTL(fnum))//'.h5'
CALL LOAD_GPE_FIELD(filename,icpoint,time,dt_argle)
! Construct psi from psire and psiim.
DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx
cgp_phys_tmp1(i1,i2_loc) = CMPLX(psire(i1,i2_loc),&
                          psiim(i1,i2_loc),KIND=GP)
ENDDO
ENDDO
#endif

END SUBROUTINE read_field


SUBROUTINE restart_from_argle_data(nrestart)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi

IMPLICIT NONE

INTEGER:: i1,i2,i3,i4
INTEGER:: i1_loc,i2_loc
INTEGER:: ishell
INTEGER:: k1,k2,ksqr,mshl
!INTEGER:: time,icpoint
INTEGER, INTENT(IN):: nrestart
INTEGER:: icpoint
INTEGER:: jtime
REAL(KIND=GP):: time
REAL(KIND=GP):: dt_argle
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky
REAL(KIND=GP):: rk,rk2
REAL(KIND=GP):: atemp,amin,amax,tmax,tmin
REAL(KIND=GP):: u1,u2,u3,u4,u5,u6,ss1
COMPLEX(KIND=GP):: cwf
!! Regrouping of data files.
CHARACTER(100)::fnn,prcjj
INTEGER:: nregrpfile,ii,jj,kk,nplnperproc_orignal

CHARACTER(LEN=100):: filename

!CALL fft(plan_vel1%pname,inverse)
!CALL fft(plan_vel2%pname,inverse)
!CALL fft(plan_vel3%pname,inverse)

WRITE(fnum,'(i8)')nrestart

#ifdef IOEVERYPROC

IF (nsplit_orignal .NE. nsplit_current .AND. nsplit_orignal .GT. nsplit_current) THEN

nregrpfile=nsplit_orignal/nsplit_current
nplnperproc_orignal = Ny/nsplit_orignal

i2_loc=0

DO ii = 1,nregrpfile
jj = nregrpfile*myrank + ii

WRITE(prcjj,'(i8)') jj

OPEN(UNIT=11,FILE='restart/wf/wf'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) icpoint,time,dt_argle

DO kk = 1,nplnperproc_orignal
i2_loc = i2_loc + 1

READ(11) (psi(i1,i2_loc),i1=1,Nx)

ENDDO

CLOSE(11)

ENDDO

ELSEIF (nsplit_orignal .NE. nsplit_current .AND. nsplit_orignal .LT. nsplit_current) THEN
 
jj = ((Ny/nsplit_current)*myrank)/(Ny/nsplit_orignal)+1
WRITE(prcjj,'(i8)') jj

OPEN(UNIT=11,FILE='restart/wf/wf'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) icpoint,time,dt_argle

i2_loc=0
DO ii = (jj-1)*(Ny/nsplit_orignal)+1, jj*(Ny/nsplit_orignal)

IF (ii .GT. (Ny/nsplit_current)*myrank .AND. ii .LE. (Ny/nsplit_current)*(myrank+1) ) THEN
i2_loc = i2_loc + 1

READ(11) (psi(i1,i2_loc),i1=1,Nx)

ELSE

READ(11) (cwf,i1=1,Nx)

ENDIF
ENDDO

CLOSE(11)

ELSE

WRITE(prc,'(i8)') (myrank+1)

OPEN(UNIT=11,FILE='restart/wf/wf'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) icpoint,time,dt_argle

DO i2_loc=1,local_Ny_cmpx
READ(11) (psi(i1,i2_loc),i1=1,Nx)
ENDDO
CLOSE(11)

ENDIF

#endif

#ifdef IOHDF5
filename = 'restart/wf/wf_'//TRIM(ADJUSTL(fnum))//'.h5'
CALL LOAD_GPE_FIELD(filename,icpoint,time,dt_argle)
! Construct psi from psire and psiim.
DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx
psi(i1,i2_loc) = CMPLX(psire(i1,i2_loc),&
                          psiim(i1,i2_loc),KIND=GP)
ENDDO
ENDDO
#endif

!tstart = jtime
timestart = time
CALL fft(plan_psi%pname,forward)

IF (myrank .EQ. 0) THEN
PRINT*,'Restart file (ARGLE):',nrestart
PRINT*,'Restart no. as inside file (ARGLE):',icpoint
PRINT*,'Restart time (ARGLE):',time
ENDIF

!! Set the time and file counter to zero.
timestart = zero
nrst = 0

END SUBROUTINE restart_from_argle_data

#endif
