!!=============================================
#ifdef SOL_GPE

SUBROUTINE diagnostics_gpe(spoint,jtime,time)

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD
USE mod_fft
USE mod_mpi
#ifdef GPPART
USE mod_gpparticles
#endif

IMPLICIT NONE

INTEGER, INTENT(IN):: spoint,jtime
REAL(KIND=GP), INTENT(IN):: time
INTEGER:: idmn
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2, ksqr, mshl
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky, rk,rk2,kradius2
REAL(KIND=GP):: enorm,e_mode,enst_mode,heli_mode,lsc_mode
REAL(KIND=GP):: ekint_mode,ekini_mode,ekinc_mode,equnt_mode
REAL(KIND=GP):: eintr_mode,eoccp_mode
REAL(KIND=GP):: egrad,ekint,equnt,eintr,echem,etot,etemp
REAL(KIND=GP):: etemp1,etemp2,etemp3,etemp4,etemp5,etemp6
REAL(KIND=GP):: momx,momy,momz,tempmom
!REAL(KIND=GP), DIMENSION(0:nshell,1:2):: spec!! Consider making this allocatable.
REAL(KIND=GP):: sptemp1_gp,sptemp2_gp,sptemp3_gp
REAL(KIND=GP):: sptemp4_gp,sptemp5_gp,sptemp6_gp
REAL(KIND=GP):: rholoc
COMPLEX(KIND=GP):: aa,aa1,aa2,aa3,bb,bb1,bb2,bb3
COMPLEX(KIND=GP):: cc,cc1,cc2,cc3,dd,dd1,dd2,dd3
COMPLEX(KIND=GP):: ee,ee1,ee2,ee3,ff,ff1,ff2,ff3
COMPLEX(KIND=GP):: gg,gg1,gg2,gg3
REAL(KIND=GP):: hh,hh1,hh2,hh3
REAL(KIND=GP):: norm, normtemp
REAL(KIND=GP), PARAMETER:: smlnum = 1e-30
#ifdef GPPART
INTEGER:: ipart
REAL(KIND=GP):: eppot,eobj_kin
REAL(KIND=GP):: gp_eobjrep_factor,erep,esrep
REAL(KIND=GP):: etemp7, etemp8
REAL(KIND=GP):: mom_obj(1:ndim)
#endif

WRITE(fnrst,'(i8)')nrst

!! Gradient of psi in Fourier space.
DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

cgp_four_tmp1(i2,i1_loc) = zi*kx*kpsi(i2,i1_loc)
cgp_four_tmp2(i2,i1_loc) = zi*ky*kpsi(i2,i1_loc)

ENDDO
ENDDO

!! Get the gradients in physical space.
CALL fft(plan_gptmp1%pname,inverse)
CALL fft(plan_gptmp2%pname,inverse)

!! Wave function in real/physical space.
!CALL fft(plan_psi%pname,inverse)

#ifdef GPPART
!! Set ipart = 0, to get the sum over all the particles.
ipart = 0
CALL translate_particle_potential(ipart)
#endif

!! Momentum 
momx = zero; momy = zero;
equnt = zero
ekint = zero

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

!momx = momx + aimag(conjg(psi(i1,i2_loc))*cgp_phys_tmp1(i1,i2_loc))
!momy = momy + aimag(conjg(psi(i1,i2_loc))*cgp_phys_tmp2(i1,i2_loc))

rholoc = ABS(psi(i1,i2_loc))**2+smlnum

aa1 = conjg(psi(i1,i2_loc))*cgp_phys_tmp1(i1,i2_loc)
aa2 = conjg(psi(i1,i2_loc))*cgp_phys_tmp2(i1,i2_loc)

momx = momx + aimag(aa1)
momy = momy + aimag(aa2)

cgp_phys_nlin(i1,i2_loc) = abs(psi(i1,i2_loc))**2

! Total kinetic energy
ekint = ekint + half*(two*alpha)**2*(AIMAG(aa1)**2+AIMAG(aa2)**2)/rholoc
                                

! Quantum pressure energy.

equnt = equnt + two*alpha**2*(REAL(aa1)**2+REAL(aa2)**2)/rholoc

#ifdef GPPART
phys_cmpxtmp(i1,i2_loc) = cgp_phys_trans_vobj(i1,i2_loc)*psi(i1,i2_loc)
#endif

ENDDO
ENDDO

momx = two*alpha*momx*dvol/(lengthx*lengthy)
momy = two*alpha*momy*dvol/(lengthx*lengthy)

!! Collect components of momentum from different procs.
CALL mpi_reduce(momx,tempmom,1,mpi_double_precision, mpi_sum,0, &
                mpi_comm_world, ierr)
IF (myrank .EQ. 0) momx = tempmom
CALL mpi_reduce(momy,tempmom,1,mpi_double_precision, mpi_sum,0, &
                mpi_comm_world, ierr)
IF (myrank .EQ. 0) momy = tempmom

!! Collect Quantum pressure energy
CALL MPI_REDUCE(equnt,etemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) equnt = etemp*dvol/(lengthx*lengthy)

!! Collect Total kineticl energy
CALL MPI_REDUCE(ekint,etemp6,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) ekint = etemp6*dvol/(lengthx*lengthy)

!! |psi|^2 in fourier space.
CALL fft(plan_gpnlin%pname,forward)

#ifdef GPPART
!! Translated psi*V(r-q) in fourier space.
CALL fft(plan_cmpxtmp%pname,forward)
#endif


!---Dealias the |psi|^2 term in Fourier-space here.
DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)

ksqr = k1**2+k2**2

IF (ksqr .GE. kalias**2) THEN

cgp_four_nlin(i2,i1_loc) = czero

#ifdef GPPART
four_cmpxtmp(i2,i1_loc) = czero
#endif

ENDIF

!! Store it for interaction energy spectra.
cgp_four_str1(i2,i1_loc) = cgp_four_nlin(i2,i1_loc)

ENDDO
ENDDO

!!---Obtain the P_G[|psi|^2] in physical space.
CALL fft(plan_gpnlin%pname,inverse)

#ifdef GPPART
CALL fft(plan_cmpxtmp%pname,inverse)
#endif

!! Energy
egrad = zero
eintr = zero
echem = zero
etot = zero
#ifdef GPPART
eppot = zero
#endif

norm = zero

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

egrad = egrad + alpha*(ABS(cgp_phys_tmp1(i1,i2_loc))**2 &
              +        ABS(cgp_phys_tmp2(i1,i2_loc))**2)

eintr = eintr + half*gstrength*REAL(cgp_phys_nlin(i1,i2_loc))**2

echem = echem + omega_chem*REAL(psi(i1,i2_loc))**2

#ifdef GPPART
eppot = eppot + REAL(CONJG(psi(i1,i2_loc))*phys_cmpxtmp(i1,i2_loc))
#endif

norm = norm + ABS(psi(i1,i2_loc))**2

ENDDO
ENDDO

!! Collect from different procs. 

CALL MPI_REDUCE(egrad,etemp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) egrad = two*alpha*etemp1*dvol/(lengthx*lengthy)

CALL MPI_REDUCE(eintr,etemp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) eintr = two*alpha*etemp2*dvol/(lengthx*lengthy)

CALL MPI_REDUCE(echem,etemp3,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) echem = two*alpha*etemp3*dvol/(lengthx*lengthy)

CALL MPI_REDUCE(etot,etemp4,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) etot = egrad + eintr !+ echem 

#ifdef GPPART
CALL MPI_REDUCE(eppot,etemp7,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) eppot = two*alpha*etemp7*dvol/(lengthx*lengthy)
#endif

CALL MPI_REDUCE(norm,normtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) norm = normtemp*dvol

#ifdef GPPART
!! Sum of particle energies.
gp_eobjrep_factor = 64.0_GP*gpdeltaE*pi**(12)*gpr0obj**(12)/lengthx**(12) 
eobj_kin = zero
esrep = zero

mom_obj(1:ndim) = zero

DO i1 = 1, gp_N_obj

DO idmn = 1, ndim
eobj_kin = eobj_kin + gp_vel_obj(i1,idmn)**2
mom_obj(idmn) = mom_obj(idmn) + gp_mass_part*gp_vel_obj(i1,idmn)
ENDDO

DO i2 = 1, gp_N_obj

x = gp_pos_obj(i1,1)-gp_pos_obj(i2,1)
y = gp_pos_obj(i1,2)-gp_pos_obj(i2,2)

erep = gp_eobjrep_factor/(-3+COS(x)+COS(y))**6
IF (i1==i2) erep = zero
esrep = esrep + erep
ENDDO
ENDDO

eobj_kin = half*gp_mass_part*eobj_kin
esrep = half*esrep !(half is there to correct double counting)

IF ( gp_N_obj == 1)  esrep = zero
#endif

IF (myrank .EQ. 0) THEN
OPEN(UNIT=1,FILE='data/energy1.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,10) time, egrad,ekint,equnt,eintr,echem,etot
10 FORMAT(E15.6,E15.6,E15.6,E15.6,E15.6,E15.6,E15.6)
CLOSE(1)

OPEN(UNIT=1,FILE='data/norm.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,20) time,norm
20 FORMAT(E15.6,E15.6)
CLOSE(1)

OPEN(UNIT=1,FILE='data/mom.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,30) time,momx,momy
30 FORMAT(E15.6,E15.6,E15.6)
CLOSE(1)

ENDIF

!!--rho^(1/2)u= 2 alpha Im (psi* grad psi )/abs(psi)

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

rholoc = ABS(psi(i1,i2_loc))**2+smlnum

aa1 = ABS(psi(i1,i2_loc))*AIMAG(CONJG(psi(i1,i2_loc))*cgp_phys_tmp1(i1,i2_loc))
aa2 = ABS(psi(i1,i2_loc))*AIMAG(CONJG(psi(i1,i2_loc))*cgp_phys_tmp2(i1,i2_loc))

cgp_phys_tmp1(i1,i2_loc) =  two*alpha*aa1/rholoc
cgp_phys_tmp2(i1,i2_loc) =  two*alpha*aa2/rholoc

!cgp_phys_tmp1(i1,i2_loc) =  two*alpha*AIMAG(CONJG(psi(i1,i2_loc))*cgp_phys_tmp1(i1,i2_loc))/ABS(psi(i1,i2_loc))
!cgp_phys_tmp2(i1,i2_loc) =  two*alpha*AIMAG(CONJG(psi(i1,i2_loc))*cgp_phys_tmp2(i1,i2_loc))/ABS(psi(i1,i2_loc))

!!--For quantum pressure term.
phys_cmpxtmp(i1,i2_loc) = ABS(psi(i1,i2_loc))

ENDDO
ENDDO

!! Get it in Fourier space.
CALL fft(plan_gptmp1%pname,forward)
CALL fft(plan_gptmp2%pname,forward)

!! Fourier transform of abs(psi)
CALL fft(plan_cmpxtmp%pname,forward)


!enorm = (one/REAL(Nx*Ny*Nz,KIND=GP))**2
enorm = (one/(REAL(Nx,KIND=GP)*REAL(Ny,KIND=GP)))**2

spec_gp = zero

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

aa1 = cgp_four_tmp1(i2,i1_loc)
aa2 = cgp_four_tmp2(i2,i1_loc)

IF (k1 .EQ. 0 .AND. k2 .EQ. 0) THEN
bb1 = czero
bb2 = czero
ELSE
bb1 = -zi*(kx/kradius2)*aa1
bb2 = -zi*(ky/kradius2)*aa2
ENDIF

bb = bb1 + bb2
!!--Compressible part
cc1 = zi*kx*bb
cc2 = zi*ky*bb
!!--Incompressible part
dd1 = aa1 - cc1
dd2 = aa2 - cc2

!!--Grad abs(psi)
ee1 = zi*kx*four_cmpxtmp(i2,i1_loc)
ee2 = zi*ky*four_cmpxtmp(i2,i1_loc)

!!--Interaction term
ff = cgp_four_str1(i2,i1_loc)

!!--Occupation spectra term
gg = kpsi(i2,i1_loc)

ekint_mode = enorm*(ABS(aa1)**2+ABS(aa2)**2)
ekinc_mode = enorm*(ABS(cc1)**2+ABS(cc2)**2)
ekini_mode = enorm*(ABS(dd1)**2+ABS(dd2)**2)
equnt_mode = enorm*(ABS(ee1)**2+ABS(ee2)**2)
eintr_mode = enorm*ABS(ff)**2
eoccp_mode = enorm*ABS(gg)**2

IF (ksqr .GE. kalias**2) THEN
ekint_mode = zero
ekini_mode = zero
ekinc_mode = zero
equnt_mode = zero
eintr_mode = zero
eoccp_mode = zero
ENDIF

spec_gp(mshl,1) = spec_gp(mshl,1) + ekint_mode*half
!spec_gp(mshl,1) = spec_gp(mshl,1) + ekint_mode*two*alpha**2
spec_gp(mshl,2) = spec_gp(mshl,2) + ekini_mode*half
spec_gp(mshl,3) = spec_gp(mshl,3) + ekinc_mode*half
spec_gp(mshl,4) = spec_gp(mshl,4) + equnt_mode*two*alpha**2
spec_gp(mshl,5) = spec_gp(mshl,5) + eintr_mode*half*gstrength*two*alpha
spec_gp(mshl,6) = spec_gp(mshl,6) + eoccp_mode

ENDDO
ENDDO

DO mshl = 0,nshell
CALL MPI_REDUCE(spec_gp(mshl,1),sptemp1_gp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(spec_gp(mshl,2),sptemp2_gp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(spec_gp(mshl,3),sptemp3_gp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(spec_gp(mshl,4),sptemp4_gp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(spec_gp(mshl,5),sptemp5_gp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
CALL MPI_REDUCE(spec_gp(mshl,6),sptemp6_gp,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
!--------
IF (myrank .EQ. 0) THEN
spec_gp(mshl,1) = sptemp1_gp
spec_gp(mshl,2) = sptemp2_gp
spec_gp(mshl,3) = sptemp3_gp
spec_gp(mshl,4) = sptemp4_gp
spec_gp(mshl,5) = sptemp5_gp
spec_gp(mshl,6) = sptemp6_gp
ENDIF
ENDDO

IF (myrank .EQ. 0) THEN
#ifdef FIXDT
IF (MOD(jtime,gsteps)==0) CALL quantglobal_gpe(jtime,time)
IF (MOD(jtime,ssteps)==0) CALL spectrum_gpe(spoint,time)
#endif
#ifdef ADAPT
IF (MOD(time,gtime) .LT. dt) CALL quantglobal_gpe(jtime,time)
IF (MOD(time,stime) .LT. dt) CALL spectrum_gpe(spoint,time)
#endif


ENDIF


!! Temporary quick check for quantum pressure energy

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

!!--For quantum pressure term.
phys_cmpxtmp(i1,i2_loc) = ABS(psi(i1,i2_loc))

ENDDO
ENDDO

!! Fourier transform of abs(psi)
CALL fft(plan_cmpxtmp%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

IF (ksqr .GE. kalias**2) THEN
four_cmpxtmp(i2,i1_loc) = czero
cgp_four_tmp1(i2,i1_loc) = czero
cgp_four_tmp2(i2,i1_loc) = czero
ENDIF

cgp_four_tmp1(i2,i1_loc) = zi*kx*four_cmpxtmp(i2,i1_loc)
cgp_four_tmp2(i2,i1_loc) = zi*ky*four_cmpxtmp(i2,i1_loc)

ENDDO
ENDDO

!! Get the gradients in physical space.
CALL fft(plan_gptmp1%pname,inverse)
CALL fft(plan_gptmp2%pname,inverse)

equnt = zero

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

equnt = equnt + two*alpha**2*(ABS(cgp_phys_tmp1(i1,i2_loc))**2 &
              +               ABS(cgp_phys_tmp2(i1,i2_loc))**2)
ENDDO
ENDDO

CALL MPI_REDUCE(equnt,etemp5,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        0,MPI_COMM_WORLD,ierr)
IF (myrank .EQ. 0) equnt = etemp5*dvol/(lengthx*lengthy)

IF (myrank .EQ. 0) THEN
OPEN(UNIT=1,FILE='data/equnt.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,40) time, equnt
40 FORMAT(E15.6,E15.6)
CLOSE(1)
ENDIF

!----

#ifdef GPPART

DO i1 = 1, gp_N_obj
CALL translate_particle_potential(i1)
DO idmn = 1, ndim
gp_frc_obj(i1,idmn) = gpfrcpart(idmn)
ENDDO
ENDDO

IF (myrank .EQ. 0) THEN
!ipart=1
DO i1 = 1, gp_N_obj

WRITE(nobj,'(i8)')i1

OPEN(UNIT=1,FILE='data/gpobj_pos_'//trim(adjustl(nobj))//'.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,45) time, gp_pos_obj(i1,1),gp_pos_obj(i1,2)
45 FORMAT(E15.6,E15.6,E15.6)
CLOSE(1)

OPEN(UNIT=1,FILE='data/gpobj_vel_'//trim(adjustl(nobj))//'.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,45) time, gp_vel_obj(i1,1),gp_vel_obj(i1,2)
46 FORMAT(E15.6,E15.6,E15.6)
CLOSE(1)

OPEN(UNIT=1,FILE='data/gpobj_frc_'//trim(adjustl(nobj))//'.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,45) time, gp_frc_obj(i1,1),gp_frc_obj(i1,2)
47 FORMAT(E15.6,E15.6,E15.6)
CLOSE(1)
ENDDO

! Write the energy.
OPEN(UNIT=1,FILE='data/egpobj.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,48) time, eppot, eobj_kin, esrep
CLOSE(1)
48 FORMAT(E15.6,E15.6,E15.6,E15.6)

OPEN(UNIT=1,FILE='data/momgpobj.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,49) time, mom_obj(1), mom_obj(2)
CLOSE(1)
49 FORMAT(E15.6,E15.6,E15.6)

ENDIF
#endif


END SUBROUTINE diagnostics_gpe

SUBROUTINE quantglobal_gpe(jtime,time)

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD

IMPLICIT NONE

INTEGER, INTENT(IN):: jtime
REAL(KIND=GP), INTENT(IN):: time
REAL(KIND=GP):: ener,enst,heli,Mlsc
REAL(KIND=GP):: ekint,ekini,ekinc,equnt,eintr
!REAL(KIND=GP), INTENT(IN),DIMENSION(0:nshell,1:2):: spec

ekint = sum(spec_gp(1:nshell,1))
ekini = sum(spec_gp(1:nshell,2))
ekinc = sum(spec_gp(1:nshell,3))
equnt = sum(spec_gp(1:nshell,4))
eintr = sum(spec_gp(0:nshell,5))

WRITE(fnrst,'(i8)')nrst

OPEN(UNIT=1,FILE='data/energy.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,10) time, ekint,ekini,ekinc,equnt,eintr
10 FORMAT(E15.6,E15.6,E15.6,E15.6,E15.6,E15.6)
CLOSE(1)

OPEN(UNIT=1,FILE='data/n0.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,20) time, spec_gp(0,6)
20 FORMAT(E15.6,E15.6)
CLOSE(1)

END SUBROUTINE quantglobal_gpe

SUBROUTINE spectrum_gpe(spoint,time)

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD

IMPLICIT NONE

INTEGER, INTENT(IN):: spoint
REAL(KIND=GP), INTENT(IN):: time
INTEGER:: i1
!REAL(KIND=GP), INTENT(IN),DIMENSION(0:nshell,1:2):: spec

WRITE(fnrst,'(i8)')nrst
WRITE(fnum,'(i8)') spoint

OPEN(UNIT=1,FILE='spectra/spec'//trim(adjustl(fnum))//'.dat',&
        STATUS='UNKNOWN',IOSTAT=ios)
DO i1=1,nshell
WRITE(1,10) dfloat(i1),spec_gp(i1,1),spec_gp(i1,2),spec_gp(i1,3),&
                       spec_gp(i1,4),spec_gp(i1,5),spec_gp(i1,6)
10 FORMAT(E15.6,E15.6,E15.6,E15.6,E15.6,E15.6,E15.6)
ENDDO
CLOSE(1)

END SUBROUTINE spectrum_gpe

#ifdef GPSTS
SUBROUTINE sts_gpe(stspoint,time)

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD
USE mod_fft
USE mod_mpi

IMPLICIT NONE

INTEGER, INTENT(IN):: stspoint
REAL(KIND=GP), INTENT(IN):: time
INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2, ksqr, mshl
REAL(KIND=GP):: kx,ky, rk,rk2,kradius2

!----

WRITE(fnrst,'(i8)')nrst
WRITE(prc,'(i8)') (myrank+1)
WRITE(fnum,'(i8)') stspoint

!IF (myrank .EQ. 0) THEN

!OPEN(UNIT=1,FILE='data/sts_time.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
!WRITE(1,10) time
!10 FORMAT(E15.6)
!CLOSE(1)

!!--(kx,0,0)
!OPEN(UNIT=1,FILE='data/sts1d/kpsi_y0z0_nrst'//TRIM(ADJUSTL(fnrst))//'.dat',&
!FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
!WRITE(1) (kpsi(i1,1,1),i1=1,Nx)
!!PRINT*, (abs(kpsi(i1,1,1)),i1=1,Nx)
!CLOSE(1)
!!--(0,kz,0)
!OPEN(UNIT=1,FILE='data/sts1d/kpsi_x0y0_nrst'//TRIM(ADJUSTL(fnrst))//'.dat',&
!FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
!WRITE(1) (kpsi(1,i3,1),i3=1,Nz)
!CLOSE(1)
!!--(kx,ky,0) --- Averaged over y-direction
!OPEN(UNIT=1,FILE='data/sts2d/kpsi_y0_nrst'//TRIM(ADJUSTL(fnrst))//'.dat',&
!FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
!WRITE(1) ((kpsi(i1,i3,1),i1=1,Nx),i3=1,Nz)
!CLOSE(1)

!ENDIF

!CALL mpi_barrier(MPI_COMM_WORLD,ierr)

!!--(kx,0,kz) --- Averaged over z-direction
!OPEN(UNIT=1,FILE='data/sts2d/kpsi_z0_nrst'//TRIM(ADJUSTL(fnrst))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
!FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
!WRITE(1) ((kpsi(i1,1,i2_loc),i1=1,Nx),i2_loc=1,local_Ny_cmpx)
!CLOSE(1)

!!--(0,ky,kz) --- Averaged over x-direction
!OPEN(UNIT=1,FILE='data/sts2d/kpsi_x0_nrst'//TRIM(ADJUSTL(fnrst))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
!FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
!WRITE(1) ((kpsi(1,i3,i2_loc),i3=1,Nz),i2_loc=1,local_Ny_cmpx)
!CLOSE(1)

!CALL mpi_barrier(MPI_COMM_WORLD,ierr)


!! With many MPI_barriers to check data writing.

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

IF (myrank .EQ. 0) THEN
OPEN(UNIT=1,FILE='data/sts_time.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
WRITE(1,10) time
10 FORMAT(E15.6)
CLOSE(1)
ENDIF

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

IF (myrank .EQ. 0) THEN
!!--(kx,0,0)
OPEN(UNIT=1,FILE='data/sts1d/kpsi_y0_nrst'//TRIM(ADJUSTL(fnrst))//'.dat',&
FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
WRITE(1) (kpsi(i1,1,1),i1=1,Nx)
!PRINT*, (abs(kpsi(i1,1,1)),i1=1,Nx)
CLOSE(1)
ENDIF

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

IF (myrank .EQ. 0) THEN
!!--(0,kz,0)
OPEN(UNIT=1,FILE='data/sts1d/kpsi_x0y0_nrst'//TRIM(ADJUSTL(fnrst))//'.dat',&
FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
WRITE(1) (kpsi(1,i3,1),i3=1,Nz)
CLOSE(1)
ENDIF

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

IF (myrank .EQ. 0) THEN
!!--(kx,ky,0) --- Averaged over y-direction
OPEN(UNIT=1,FILE='data/sts2d/kpsi_y0_nrst'//TRIM(ADJUSTL(fnrst))//'.dat',&
FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
WRITE(1) ((kpsi(i1,i3,1),i1=1,Nx),i3=1,Nz)
CLOSE(1)
ENDIF

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

!!--(kx,0,kz) --- Averaged over z-direction
OPEN(UNIT=1,FILE='data/sts2d/kpsi_z0_nrst'//TRIM(ADJUSTL(fnrst))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
WRITE(1) ((kpsi(i1,1,i2_loc),i1=1,Nx),i2_loc=1,local_Ny_cmpx)
CLOSE(1)

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

!!--(0,ky,kz) --- Averaged over x-direction
OPEN(UNIT=1,FILE='data/sts2d/kpsi_x0_nrst'//TRIM(ADJUSTL(fnrst))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED', status='UNKNOWN', POSITION='APPEND', IOSTAT=ios)
WRITE(1) ((kpsi(1,i3,i2_loc),i3=1,Nz),i2_loc=1,local_Ny_cmpx)
CLOSE(1)

CALL mpi_barrier(MPI_COMM_WORLD,ierr)

END SUBROUTINE sts_gpe
#endif

#endif
