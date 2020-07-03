#ifdef SOL_GPE

#ifdef GPE_ARGLE

SUBROUTINE ARGLE(time)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
USE mtmod
#ifdef GPPART
USE mod_gpparticles
#endif

IMPLICIT NONE
REAL(KIND=GP), INTENT(IN):: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: Nsquare
INTEGER:: mshl
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius2
REAL(KIND=GP):: uadx,uady
COMPLEX(KIND=GP):: prefacrhs
COMPLEX(KIND=GP):: aa1,aa2
COMPLEX(KIND=GP):: tmpstr1,tmpstr2
COMPLEX(KIND=GP):: noise_term
REAL(KIND=GP):: ran1,ran2,v1,v2,rsq
#ifdef GPPART
INTEGER:: ipart
COMPLEX(KIND=GP):: potterm
#endif

!prefacrhs = (one-zi*disp_factor)/zi

!! Updated wave function in real/physical space.
CALL fft(plan_psi%pname,inverse)

!! Note: It is best to compute the contribution of the
!! particles to the GP Equation. Everything is contained
!! in the V(x-q), i.e. potential term. This information
!! will be held in cgp_four_trans_vobj.
!! This will also avoid any possible mix up in the use
!! of temporary arrays, usually used to calculate gradients.
#ifdef GPPART
ipart = 0
CALL translate_particle_potential(ipart)
#endif

!! |psi|^2 in physical space.
DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx
cgp_phys_nlin(i1,i2_loc) = abs(psi(i1,i2_loc))**2
ENDDO
ENDDO

!! |psi|^2 in fourier space.
CALL fft(plan_gpnlin%pname,forward)


!! Linear terms and Dealias the |psi|^2 term in Fourier-space.
DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

ksqr = k1**2+k2**2
kradius2 = kx**2 + ky**2

#if defined(GPE_ADVECTION) || defined(GPE_COUNTERFLOW)
cgp_four_tmp1(i2,i1_loc) = zi*kx*kpsi(i2,i1_loc)
cgp_four_tmp2(i2,i1_loc) = zi*ky*kpsi(i2,i1_loc)
#endif

!Dealias the |psi|^2 term
IF (ksqr .GE. kalias**2) THEN
cgp_four_nlin(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

#if defined(GPE_ADVECTION) || defined(GPE_COUNTERFLOW)
!! Get the gradients in physical space.
CALL fft(plan_gptmp1%pname,inverse)
CALL fft(plan_gptmp2%pname,inverse)
#endif

!!---Obtain the P_G[|psi|^2] in physical space.
CALL fft(plan_gpnlin%pname,inverse)

!!---Now calculate [P_G|psi|^2]psi term
DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
y = (i2-1)*dy
DO i1 = 1, Nx
x = (i1-1)*dx

cgp_phys_nlin(i1,i2_loc) = -gstrength*cgp_phys_nlin(i1,i2_loc)*psi(i1,i2_loc)

#ifdef GPPART
potterm = -cgp_phys_trans_vobj(i1,i2_loc)*psi(i1,i2_loc)
cgp_phys_nlin(i1,i2_loc) = cgp_phys_nlin(i1,i2_loc) + potterm
#endif

#ifdef GPE_ADVECTION
! Advective part.
uadx =  SIN(x)*COS(y)
uady = -COS(x)*SIN(y)

aa1 = cgp_phys_tmp1(i1,i2_loc)
aa2 = cgp_phys_tmp2(i1,i2_loc)

!cgp_phys_tmp1(i1,i2_loc) = czero
!cgp_phys_tmp2(i1,i2_loc) = czero

!!cgp_phys_tmp1(i1,i2,i3_loc) = -zi*(uadx*aa1+uady*aa2+uadz*aa3)
!!cgp_phys_tmp2(i1,i2,i3_loc) = -(half**2)*(one/alpha)*(uadx**2+uady**2+uadz**2)*psi(i1,i2,i3_loc)
!! Save FFT operation by combining these.
cgp_phys_nlin(i1,i2_loc) = cgp_phys_nlin(i1,i2_loc) &
        -zi*(uadx*aa1+uady*aa2) &
        -(half**2)*(one/alpha)*(uadx**2+uady**2)*psi(i1,i2_loc)

#endif 
! Above is the end of advection.

#ifdef GPE_COUNTERFLOW
! Constant counterflow velocity.
uadx = vnx
uady = vny

aa1 = cgp_phys_tmp1(i1,i2_loc)
aa2 = cgp_phys_tmp2(i1,i2_loc)

cgp_phys_nlin(i1,i2_loc) = cgp_phys_nlin(i1,i2_loc) &
        -zi*(uadx*aa1+uady*aa2)

#endif


#ifdef GPE_SGLE
!! SGLE Noise term.
DO 
!ran1 = grnd(); ran2 = grnd()
CALL RANDOM_NUMBER(ran1)
CALL RANDOM_NUMBER(ran2)
v1 = two*ran1 - one
v2 = two*ran2 - one
rsq = v1**2 + v2**2
IF (rsq .GT. zero .AND. rsq .LT. one) EXIT
ENDDO
rsq = SQRT(-two*LOG(rsq)/rsq)
v1 = v1*rsq
v2 = v2*rsq
phys_cmpxtmp(i1,i2_loc) = CMPLX(v1,v2)*(one/SQRT(dvol))*noisesigma
#endif

ENDDO
ENDDO

!!---Obtain the [P_G|psi|^2]psi+terms in Fourier space
CALL fft(plan_gpnlin%pname,forward)

!!---For Stochastic Ginzburg Landau equation (SGLE)
!! Get the chemical potential from the adhoc equation.
#ifdef GPE_SGLE
CALL finiteT_adhoc_chempot(time,omega_chem)
!! Obtain the noise in Fourier space.
CALL fft(plan_cmpxtmp%pname,forward)
#endif


!! Update kpsi
DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

ksqr = k1**2+k2**2
kradius2 = kx**2 + ky**2

#ifndef GPE_SGLE

tmpstr1 = cgp_four_nlin(i2,i1_loc) + omega_chem*kpsi(i2,i1_loc)
tmpstr2 = (one+alpha*kradius2*dt)
kpsi(i2,i1_loc) = (kpsi(i2,i1_loc)+dt*tmpstr1)/tmpstr2

#else

noise_term = four_cmpxtmp(i2,i1_loc)*SQRT(dt)
!! Dealias it.
IF (ksqr .GE. kalias**2) THEN
noise_term = czero
ENDIF

!tmpstr1 = cgp_four_nlin(i2,i1_loc)
!tmpstr2 = one + (alpha*kradius2 - omega_chem)*dt
tmpstr1 = cgp_four_nlin(i2,i1_loc) + omega_chem*kpsi(i2,i1_loc)
tmpstr2 = one + (alpha*kradius2 )*dt
kpsi(i2,i1_loc) = (kpsi(i2,i1_loc)+dt*tmpstr1 + noise_term)/tmpstr2 !+ noise_term

#ifdef GPE_SELFTRUNC
IF (ksqr .GT. ktrunc**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF
#endif

#endif

ENDDO
ENDDO

END SUBROUTINE ARGLE

SUBROUTINE finiteT_adhoc_chempot(time,chempot)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
!
IMPLICIT NONE

REAL(KIND=GP), INTENT(IN):: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: Nsqaure
INTEGER:: mshl,nn
REAL(KIND=GP):: x,y
REAL(KIND=GP):: norm,tempnorm
REAL(KIND=GP), INTENT(INOUT):: chempot

! Compute the current norm.
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

chempot = chempot - dt*(nuN/(lengthx*lengthy))*(norm-Nparticle)

!!IF (myrank .EQ. 0) write(24,*) chempot
!IF (myrank .EQ. 0) THEN
!OPEN(UNIT=1,FILE='data/cpot_sgle.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
!WRITE(1,10) time,chempot
!10 FORMAT(E15.6,E15.6)
!CLOSE(1)
!ENDIF

END SUBROUTINE finiteT_adhoc_chempot

SUBROUTINE write_chempot_for_restart(time,chempot)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
!
IMPLICIT NONE

REAL(KIND=GP), INTENT(IN):: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: Nsqaure
INTEGER:: mshl,nn
REAL(KIND=GP):: x,y
REAL(KIND=GP):: norm,tempnorm
REAL(KIND=GP), INTENT(IN):: chempot

WRITE(fnrst,'(i8)')nrst
IF (myrank .EQ. 0) THEN
OPEN(UNIT=1,FILE='data/cpot_sgle.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,5) time,chempot
5 FORMAT(E15.6,E15.6)
CLOSE(1)
ENDIF

END SUBROUTINE write_chempot_for_restart

#else

SUBROUTINE gpe_rk4

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
#ifdef GPPART
USE mod_gpparticles
#endif

IMPLICIT NONE
INTEGER:: i1,i2,i3
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,k3,ksqr
INTEGER:: Nsquare 
INTEGER:: mshl
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius2
#ifdef GPPART
INTEGER:: rkpartstage
#endif

!! Assumes: kpsi contains fourier transform of psi.

!! Store the FT of psi.
DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

cgp_four_str1(i2,i1_loc) = kpsi(i2,i1_loc)
! Important to set to zero.
cgp_four_str2(i2,i1_loc) = czero

ENDDO
ENDDO

! k1
CALL gpe_rhs

!---Particle
#ifdef GPPART
rkpartstage = 1
CALL gp_particle_motion_rk4(rkpartstage)
#endif 

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

!! Accumulate the RK4 terms.
cgp_four_str2(i2,i1_loc) = cgp_four_str2(i2,i1_loc) &
        + dt*oney6*cgp_four_nlin(i2,i1_loc)

!! Update the wave function.
kpsi(i2,i1_loc) = cgp_four_str1(i2,i1_loc) &
        + dt*half*cgp_four_nlin(i2,i1_loc)

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

!! Updated wave function in real/physical space.
CALL fft(plan_psi%pname,inverse)

! K2
CALL gpe_rhs

!---Particle
#ifdef GPPART
rkpartstage = 2
CALL gp_particle_motion_rk4(rkpartstage)
#endif

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

!! Accumulate the RK4 terms. Add K2
cgp_four_str2(i2,i1_loc) = cgp_four_str2(i2,i1_loc) &
        + dt*oney6*two*cgp_four_nlin(i2,i1_loc)
!! Update the wave function.
kpsi(i2,i1_loc) = cgp_four_str1(i2,i1_loc) &
        + dt*half*cgp_four_nlin(i2,i1_loc)

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

!! Updated wave function in real/physical space.
CALL fft(plan_psi%pname,inverse)

! K3
CALL gpe_rhs

!---Particle
#ifdef GPPART
rkpartstage = 3
CALL gp_particle_motion_rk4(rkpartstage)
#endif

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

!! Accumulate the RK4 terms. Add K3
cgp_four_str2(i2,i1_loc) = cgp_four_str2(i2,i1_loc) &
        + dt*oney6*two*cgp_four_nlin(i2,i1_loc)
!! Update the wave function.
kpsi(i2,i1_loc) = cgp_four_str1(i2,i1_loc) &
        + dt*cgp_four_nlin(i2,i1_loc)

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

!! Updated wave function in real/physical space.
CALL fft(plan_psi%pname,inverse)

! K4
CALL gpe_rhs

!---Particle
#ifdef GPPART
rkpartstage = 4
CALL gp_particle_motion_rk4(rkpartstage)
#endif

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

!! Accumulate the RK4 terms. Add K3
cgp_four_str2(i2,i1_loc) = cgp_four_str2(i2,i1_loc) &
        + dt*oney6*cgp_four_nlin(i2,i1_loc)
!! Final update of the wave function gives its values at the next time step.
kpsi(i2,i1_loc) = cgp_four_str1(i2,i1_loc) &
        + cgp_four_str2(i2,i1_loc)

IF (ksqr .GE. kalias**2) THEN
kpsi(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

! Update the physical space array of psi
CALL fft(plan_psi%pname,inverse)

END SUBROUTINE gpe_rk4

!!------------------------------------------------------------------------

SUBROUTINE gpe_rhs

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
USE mtmod
#ifdef GPPART
USE mod_gpparticles
#endif

IMPLICIT NONE
INTEGER:: i1,i2,i3
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
INTEGER:: Ncube
INTEGER:: mshl
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius2
COMPLEX(KIND=GP):: prefacrhs
COMPLEX(KIND=GP):: gpdispterm,gpfrcterm
REAL(KIND=GP):: ran
#ifdef GPPART
INTEGER:: ipart
COMPLEX(KIND=GP):: potterm
#endif
CHARACTER(100)::fnn,prcjj

prefacrhs = (one-zi*disp_factor)/zi

!! Updated wave function in real/physical space.
!CALL fft(plan_psi%pname,inverse)

!! Note: It is best to compute the contribution of the
!! particles to the GP Equation. Everything is contained
!! in the V(x-q), i.e. potential term. This information
!! will be held in cgp_four_trans_vobj.
!! This will also avoid any possible mix up in the use
!! of temporary arrays, usually used to calculate gradients.
#ifdef GPPART
ipart = 0
CALL translate_particle_potential(ipart)
#endif

!! |psi|^2 in physical space.
DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx
cgp_phys_nlin(i1,i2_loc) = abs(psi(i1,i2_loc))**2
ENDDO
ENDDO

!! |psi|^2 in fourier space.
CALL fft(plan_gpnlin%pname,forward)


!! Linear terms and Dealias the |psi|^2 term in Fourier-space.
DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

ksqr = k1**2+k2**2
kradius2 = kx**2 + ky**2

cgp_four_tmp1(i2,i1_loc) = (alpha*kradius2-omega_chem)*kpsi(i2,i1_loc)

!Dealias the |psi|^2 term
IF (ksqr .GE. kalias**2) THEN
cgp_four_nlin(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

!!---Obtain the P_G[|psi|^2] in physical space.
CALL fft(plan_gpnlin%pname,inverse)

!!---Now calculate g[P_G|psi|^2]psi term

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx
cgp_phys_nlin(i1,i2_loc) = gstrength*cgp_phys_nlin(i1,i2_loc)*psi(i1,i2_loc)

#ifdef GPPART
potterm = cgp_phys_trans_vobj(i1,i2_loc)*psi(i1,i2_loc)
cgp_phys_nlin(i1,i2_loc) = cgp_phys_nlin(i1,i2_loc) + potterm
!print*,'hello',ABS(cgp_phys_trans_vobj(i1,i2_loc)*psi(i1,i2_loc))
#endif

ENDDO
ENDDO

!!----Check the translation of the potential
!WRITE(fnn,'(i8)') 1
!WRITE(prcjj,'(i8)') myrank+1
!OPEN(UNIT=11,FILE='Visu/Pot'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
!        status='UNKNOWN', IOSTAT=ios)
!WRITE(11,*) ((REAL(cgp_phys_trans_vobj(i1,i2_loc)*ABS(psi(i1,i2_loc))**2),i1=1,Nx),i2_loc=1, local_Ny_cmpx)
!CLOSE(11)
!!


!!---Obtain the [P_G|psi|^2]psi in Fourier space
CALL fft(plan_gpnlin%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

ksqr = k1**2+k2**2
kradius2 = kx**2 + ky**2

cgp_four_nlin(i2,i1_loc) = prefacrhs*(cgp_four_tmp1(i2,i1_loc) &
                           + cgp_four_nlin(i2,i1_loc))

#ifdef GPE_FRCDISP

!! Friction
!IF (ksqr .LT. 3.0**2) THEN
!four_cmpxtmp(i2,i1_loc) = four_cmpxtmp(i2,i1_loc) &
!        - friclowk_gp*kpsi(i2,i1_loc)
!ENDIF

!! Hypoviscosity
IF (ksqr .NE. 0) THEN
cgp_four_nlin(i2,i1_loc) = cgp_four_nlin(i2,i1_loc) &
         - vishypo_gp*kradius2**(-8)*kpsi(i2,i1_loc)
ENDIF

!! Hyperviscosity
!IF (ksqr .GE. 16.0d0) THEN
cgp_four_nlin(i2,i1_loc) = cgp_four_nlin(i2,i1_loc) &
        - vishype_gp*kradius2**8*kpsi(i2,i1_loc)
!
!ENDIF

IF (ksqr .GE. kf1_gp**2 .AND. ksqr .LE. kf2_gp**2) THEN
CALL RANDOM_NUMBER(ran)
!ran = grnd()
ran = two*pi*ran

cgp_four_nlin(i2,i1_loc) = cgp_four_nlin(i2,i1_loc) &
                    + famp_gp*EXP(zi*ran)

!four_cmpxtmp(i2,i1_loc) = four_cmpxtmp(i2,i1_loc) &
!                    +  famp_gp*EXP(zi*ran)*EXP(-((SQRT(kradius2)-k0)**2)/(two*stndev**2))
ENDIF

!CALL gpe_dissipation(i2,i1_loc,k1,k2,gpdispterm)
!CALL gpe_force(k1,k2,gpfrcterm)

!four_cmpxtmp(i2,i1_loc) = four_cmpxtmp(i2,i1_loc) &
!            + gpdispterm + gpfrcterm
#endif

!!---Dealias the RHS term in Fourier-space here.
IF (ksqr .GE. kalias**2) THEN
cgp_four_nlin(i2,i1_loc) = czero
ENDIF

ENDDO
ENDDO

END SUBROUTINE gpe_rhs

#endif
#endif
