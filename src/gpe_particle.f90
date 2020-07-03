#ifdef GPPART

SUBROUTINE particle_potential

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
USE mod_gpparticles

IMPLICIT NONE

INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius2
REAL(KIND=GP):: c1,c2,c3, sx,sy,sz, xprd,yprd,zprd
REAL(KIND=GP):: rprd
REAL(KIND=GP):: VGaussVal 

!ppotd = 2.0_GP
!ppotV0 = 100.0_GP
!ppotthick = 2.0_GP


!ppotV0 = ppotV0*gstrength!*dexp(one/4.0_dbl)
!ppotd = ppotd*xi
!ppotthick = ppotthick*xi

c1 = 2*pi/(two*lengthx)
c2 = 2*pi/(two*lengthy)

!c1 = 2*pi/(1024.0_GP*lengthx)
!c2 = 2*pi/(1024.0_GP*lengthy)

DO i2_loc = 1, local_Ny_cmpx
i2 = i2_loc + local_i2_offset_cmpx
y = (i2-1)*dy
DO i1 = 1, Nx
x = (i1-1)*dx

sx = SIN(c1*x)
sy = SIN(c2*y)

xprd = (one/c1)*(sx + oney6*sx**3 + (one/40.0_GP)*sx**5)
yprd = (one/c2)*(sy + oney6*sy**3 + (one/40.0_GP)*sy**5)

rprd = SQRT(xprd**2+yprd**2)

VGaussVal = ppotV0*EXP(-rprd**2/(two*ppotd**2))

phys_cmpxtmp(i1,i2_loc) = VGaussVal

ENDDO
ENDDO


!WRITE(334+myrank,*) ((REAL(phys_cmpxtmp(i1,i2_loc)),i1=1,Nx),i2_loc=1, local_Ny_cmpx)

!!---Obtain the particle potential in Fourier space.
CALL fft(plan_cmpxtmp%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

! where is ksqr?
! Now added below.
ksqr = k1**2+k2**2

!!---Dealias in Fourier-space here.
IF (ksqr .GE. kalias**2) THEN
four_cmpxtmp(i2,i1_loc) = czero
ENDIF

cgp_four_vobj(i2,i1_loc) = four_cmpxtmp(i2,i1_loc)
four_cmpxtmp(i2,i1_loc) = czero

ENDDO
ENDDO

END SUBROUTINE particle_potential

SUBROUTINE initialize_particles

USE mod_precision        
USE mod_constants        
USE mod_grid
USE mod_gpparticles

IMPLICIT NONE

INTEGER:: i1,i2
INTEGER:: iopoint
REAL(KIND=GP):: ptime


!gp_pos_obj(1,1) = pi
!gp_pos_obj(2,1) = pi

!gp_vel_obj(1,1) = 0.0_GP
!gp_vel_obj(2,1) = zero

!gp_f0ext_obj(1,1) = 0.5_GP!zero
!gp_f0ext_obj(2,1) = zero


!! Let every processor read the initial values.
!! ptime is discared. Starting time is taken from the field file.
WRITE(fnrst,'(i8)')nrst
! -- Position
OPEN(UNIT=1,FILE='restart/part/gppos.dat_nrst'//trim(adjustl(fnrst)))
DO i1 = 1,gp_N_obj
READ(1,*) iopoint, ptime, gp_pos_obj(1,i1), gp_pos_obj(2,i1)
!10 FORMAT(I6,E15.7,E15.7,E15.7)
ENDDO
CLOSE(1)
! -- Velocity
OPEN(UNIT=1,FILE='restart/part/gpvel.dat_nrst'//trim(adjustl(fnrst)))
DO i1 = 1,gp_N_obj
READ(1,*) iopoint, ptime, gp_vel_obj(1,i1), gp_vel_obj(2,i1)
!20 FORMAT(I6,E15.7,E15.7,E15.7)
ENDDO
CLOSE(1)
! -- Force
OPEN(UNIT=1,FILE='restart/part/gpf0ext.dat_nrst'//trim(adjustl(fnrst)))
DO i1 = 1,gp_N_obj
READ(1,*) iopoint, ptime, gp_f0ext_obj(1,i1), gp_f0ext_obj(2,i1)
!30 FORMAT(I6,E15.7,E15.7,E15.7)
ENDDO
CLOSE(1)

END SUBROUTINE initialize_particles

SUBROUTINE gppart_checkpoint(iopoint,jtime,ptime)

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_mpi
USE mod_gpparticles

IMPLICIT NONE   

INTEGER, INTENT(IN):: jtime
INTEGER, INTENT(IN):: iopoint
REAL(KIND=GP), INTENT(IN):: ptime
INTEGER:: i1

WRITE(fnrst,'(i8)')nrst

IF (myrank .EQ. 0) THEN
! -- Position
OPEN(UNIT=1,FILE='data/gppos_ckp.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
DO i1 = 1,gp_N_obj
WRITE(1,*) iopoint, ptime, gp_pos_obj(1,i1), gp_pos_obj(2,i1)
!10 FORMAT(I6,E15.7,E15.7,E15.7)
ENDDO
CLOSE(1)
! -- Velocity
OPEN(UNIT=1,FILE='data/gpvel_ckp.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
DO i1 = 1,gp_N_obj
WRITE(1,*) iopoint, ptime, gp_vel_obj(1,i1), gp_vel_obj(2,i1)
!20 FORMAT(I6,E15.7,E15.7,E15.7)
ENDDO
CLOSE(1)
! -- Force
OPEN(UNIT=1,FILE='data/gpf0ext_ckp.dat_nrst'//trim(adjustl(fnrst)),POSITION='APPEND')
DO i1 = 1,gp_N_obj
WRITE(1,*) iopoint, ptime, gp_f0ext_obj(1,i1), gp_f0ext_obj(2,i1)
!30 FORMAT(I6,E15.7,E15.7,E15.7)
ENDDO
CLOSE(1)
ENDIF

END SUBROUTINE gppart_checkpoint

SUBROUTINE translate_particle_potential(particle_trans)

!! (1) This subroutine computes the translated potential for a particular
!! particle and the force on it.
!! (2) Also computes the translated potential for all the particles when
!! needed for energy computation etc. 

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
USE mod_gpparticles

IMPLICIT NONE

INTEGER:: i1,i2
INTEGER:: i1_loc,i2_loc
INTEGER:: k1,k2,ksqr
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky,rk,rk2,kradius2
INTEGER, INTENT(IN):: particle_trans
INTEGER:: idmn,ipart
COMPLEX(KIND=GP):: sum_kq, ekq
REAL(KIND=GP):: tempgpfrcpart(1:ndim)
CHARACTER(100)::fnn,prcjj

ipart = particle_trans

!PRINT*,'particle_trans',particle_trans,'ipart',ipart,'myrank',myrank

IF (particle_trans .GE. 1 .AND. particle_trans .LE. gp_N_obj) THEN

sum_kq = czero

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

sum_kq = EXP(-zi*(kx*gp_pos_obj(1,ipart)+ky*gp_pos_obj(2,ipart)))

cgp_four_trans_vobj(i2,i1_loc) = sum_kq*cgp_four_vobj(i2,i1_loc)

ENDDO
ENDDO


!! We want trans_vk_obj in physical space.
CALL fft(plan_gptransvobj%pname,inverse)

!!----Check the translation of the potential
!WRITE(fnn,'(i8)') 1
!WRITE(prcjj,'(i8)') myrank+1
!OPEN(UNIT=11,FILE='Visu/Pot'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
!        status='UNKNOWN', IOSTAT=ios)
!WRITE(11,*) ((REAL(cgp_phys_trans_vobj(i1,i2_loc)),i1=1,Nx),i2_loc=1, local_Ny_cmpx)
!CLOSE(11)
!!-------------------------------------------

!! Calculate the shortrange repulsive forces due to 1/r^{12}.
CALL shortrange_repulsion

!! Compute |psi|^2 in physical space.
DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

phys_cmpxtmp(i1,i2_loc) = abs(psi(i1,i2_loc))**2

ENDDO
ENDDO

!! |psi|^2 in fourier space.
CALL fft(plan_cmpxtmp%pname,forward)

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

!Dealias the |psi|^2 term
IF (ksqr .GE. kalias**2) THEN
four_cmpxtmp(i2,i1_loc) = czero
ENDIF

cgp_four_tmp1(i2,i1_loc) = zi*kx*four_cmpxtmp(i2,i1_loc)
cgp_four_tmp2(i2,i1_loc) = zi*ky*four_cmpxtmp(i2,i1_loc)

ENDDO
ENDDO

!! Get the gradients in physical space.
CALL fft(plan_gptmp1%pname,inverse)
CALL fft(plan_gptmp2%pname,inverse)


!!gpfrcpart(1:ndim) = zero 
gpfrcpart(1) = zero
gpfrcpart(2) = zero

DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx

gpfrcpart(1) = gpfrcpart(1) - REAL(two*alpha*cgp_phys_tmp1(i1,i2_loc)*cgp_phys_trans_vobj(i1,i2_loc))
gpfrcpart(2) = gpfrcpart(2) - REAL(two*alpha*cgp_phys_tmp2(i1,i2_loc)*cgp_phys_trans_vobj(i1,i2_loc))

ENDDO
ENDDO

DO idmn = 1, ndim
gpfrcpart(idmn) = gpfrcpart(idmn)*dvol
ENDDO

!! Collect from all the processes

CALL MPI_REDUCE(gpfrcpart,tempgpfrcpart,ndim,mpi_double_precision,mpi_sum,&
        0,MPI_COMM_WORLD,ierr)

!! Note gp_shortrange_obj is same on all procs.
IF (myrank .EQ. 0) THEN
DO idmn = 1,ndim
gpfrcpart(idmn) = tempgpfrcpart(idmn) + gp_shortrange_obj(idmn,ipart)
gpfrcpart(idmn) = gpfrcpart(idmn) + gp_f0ext_obj(idmn,ipart)
ENDDO
ENDIF
!! Later need to the external force.
CALL MPI_BCAST(gpfrcpart,ndim,mpi_double_precision,0,MPI_COMM_WORLD,ierr)

!ENDIF

ENDIF

!!CASE (particle_trans .EQ. 0)
IF (particle_trans ==0) THEN

!if (myrank==0) print*,particle_trans, 'hello'
!PRINT*,gp_pos_obj(1,1),gp_pos_obj(2,1)

sum_kq = czero

DO i1_loc = 1, local_Nx_cmpx
i1 = i1_loc + local_i1_offset_cmpx
k1 = kwn1_c2c(i1)
kx = k1*kdx
DO i2 = 1, Ny
k2 = kwn2_c2c(i2)
ky = k2*kdy

sum_kq = czero
DO ipart = 1,gp_N_obj
sum_kq = sum_kq + EXP(-zi*(kx*gp_pos_obj(1,ipart)+ky*gp_pos_obj(2,ipart)))
!sum_kq = sum_kq + EXP(-zi*(kx*pi/two+ky*pi/two))
ENDDO

!four_cmpxtmp(i2,i1_loc) = sum_kq*cgp_four_vobj(i2,i1_loc)
cgp_four_trans_vobj(i2,i1_loc) = sum_kq*cgp_four_vobj(i2,i1_loc)
ENDDO
ENDDO


!! We want trans_vk_obj in physical space.
!CALL fft(plan_cmpxtmp%pname,inverse)
CALL fft(plan_gptransvobj%pname,inverse)

!DO i2_loc = 1, local_Ny_cmpx
!DO i1 = 1, Nx
!cgp_phys_trans_vobj(i1,i2_loc) = phys_cmpxtmp(i1,i2_loc)
!phys_cmpxtmp(i1,i2_loc) = czero
!ENDDO
!ENDDO

!!----Check the translation of the potential
!WRITE(fnn,'(i8)') 1
!WRITE(prcjj,'(i8)') myrank+1
!OPEN(UNIT=11,FILE='Visu/Pot'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
!        status='UNKNOWN', IOSTAT=ios)
!WRITE(11,*) ((REAL(cgp_phys_trans_vobj(i1,i2_loc)),i1=1,Nx),i2_loc=1, local_Ny_cmpx)
!CLOSE(11)
!!

ENDIF

END SUBROUTINE translate_particle_potential

SUBROUTINE shortrange_repulsion

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_gpparticles
USE mod_2dHD

IMPLICIT NONE

INTEGER:: i1,i2,i3
INTEGER:: i1_loc,i2_loc
REAL(KIND=GP):: x,y,z
INTEGER:: k1,k2,k3,ksqr
REAL(KIND=GP):: frepul(1:ndim),repulfactor(1:ndim),r0obj12

! Just for time being. Should set in input file and after thought.
gpdeltaE = 0.062
gpr0obj = 4.0*xi

r0obj12 = gpr0obj**(12)

repulfactor(1) = 384.0_GP*gpdeltaE*pi**(12)*r0obj12/lengthx**(12)
repulfactor(2) = 384.0_GP*gpdeltaE*pi**(12)*r0obj12/lengthx**(12)

gp_shortrange_obj(1:ndim,1:gp_N_obj) = zero
!gp_shortrange_obj(2) = zero
!gp_shortrange_obj(3) = zero

DO i1 = 1, gp_N_obj
DO i2 = 1, gp_N_obj

x = gp_pos_obj(1,i1)-gp_pos_obj(1,i2)
y = gp_pos_obj(2,i1)-gp_pos_obj(2,i2)

frepul(1) = -repulfactor(1)*SIN(x)/(-3+COS(x)+COS(y))**7
frepul(2) = -repulfactor(2)*SIN(y)/(-3+COS(x)+COS(y))**7

if (i1==i2) frepul(1:2) = zero

gp_shortrange_obj(1,i1) = gp_shortrange_obj(1,i1) + frepul(1)
gp_shortrange_obj(2,i1) = gp_shortrange_obj(2,i1) + frepul(2)

ENDDO
ENDDO

IF (gp_N_obj==1) gp_shortrange_obj(1,1) = zero
IF (gp_N_obj==1) gp_shortrange_obj(2,1) = zero

END SUBROUTINE shortrange_repulsion

SUBROUTINE gp_particle_motion_rk4(rkstage)

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_gpparticles

IMPLICIT NONE

INTEGER:: idmn
INTEGER:: ipart
INTEGER, INTENT(IN):: rkstage

SELECT CASE (rkstage)

!!-----

CASE(1)

!! At this step we have to store the current velocity and positions.
DO ipart = 1, gp_N_obj
DO idmn = 1, ndim
gp_temp_vel_obj(idmn,ipart) = gp_vel_obj(idmn,ipart)
gp_temp_pos_obj(idmn,ipart) = gp_pos_obj(idmn,ipart)
ENDDO
ENDDO

DO ipart = 1, gp_N_obj
CALL translate_particle_potential(ipart)

DO idmn = 1, ndim

k1_ugpp(idmn,ipart) = dt*gpfrcpart(idmn)/gp_mass_part
k1_qgpp(idmn,ipart) = dt*gp_vel_obj(idmn,ipart)!/gp_mass_part

gp_vel_obj(idmn,ipart) = gp_temp_vel_obj(idmn,ipart) + half*k1_ugpp(idmn,ipart)
gp_pos_obj(idmn,ipart) = gp_temp_pos_obj(idmn,ipart) + half*k1_qgpp(idmn,ipart)

ENDDO
ENDDO

!!-----
CASE(2)

DO ipart = 1, gp_N_obj
CALL translate_particle_potential(ipart)

DO idmn = 1, ndim

k2_ugpp(idmn,ipart) = dt*gpfrcpart(idmn)/gp_mass_part
k2_qgpp(idmn,ipart) = dt*gp_vel_obj(idmn,ipart)!/gp_mass_part

gp_vel_obj(idmn,ipart) = gp_temp_vel_obj(idmn,ipart) + half*k2_ugpp(idmn,ipart)
gp_pos_obj(idmn,ipart) = gp_temp_pos_obj(idmn,ipart) + half*k2_qgpp(idmn,ipart)

ENDDO
ENDDO

!!-----
CASE(3)

DO ipart = 1, gp_N_obj
CALL translate_particle_potential(ipart)

DO idmn = 1, ndim

k3_ugpp(idmn,ipart) = dt*gpfrcpart(idmn)/gp_mass_part
k3_qgpp(idmn,ipart) = dt*gp_vel_obj(idmn,ipart)!/gp_mass_part

gp_vel_obj(idmn,ipart) = gp_temp_vel_obj(idmn,ipart) + k3_ugpp(idmn,ipart)
gp_pos_obj(idmn,ipart) = gp_temp_pos_obj(idmn,ipart) + k3_qgpp(idmn,ipart)

ENDDO
ENDDO

!!-----
CASE(4)

DO ipart = 1, gp_N_obj
CALL translate_particle_potential(ipart)

DO idmn = 1, ndim

k4_ugpp(idmn,ipart) = dt*gpfrcpart(idmn)/gp_mass_part
k4_qgpp(idmn,ipart) = dt*gp_vel_obj(idmn,ipart)!/gp_mass_part


gp_vel_obj(idmn,ipart) = gp_temp_vel_obj(idmn,ipart) + oney6*&
                        (k1_ugpp(idmn,ipart) + two*k2_ugpp(idmn,ipart) &
                       + two*k3_ugpp(idmn,ipart) + k4_ugpp(idmn,ipart) ) 
    
gp_pos_obj(idmn,ipart) = gp_temp_pos_obj(idmn,ipart) + oney6*&
                        (k1_qgpp(idmn,ipart) + two*k2_qgpp(idmn,ipart) &
                       + two*k3_qgpp(idmn,ipart) + k4_qgpp(idmn,ipart) )
!!! Store the force at the 4th step.
!!gp_frc_obj(idmn,ipart) = gpfrcpart(idmn)

ENDDO
ENDDO

END SELECT

END SUBROUTINE gp_particle_motion_rk4

#endif
