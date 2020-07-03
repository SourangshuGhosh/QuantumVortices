#ifdef GPPART

MODULE mod_gpparticles

USE mod_precision
USE mod_grid
USE mod_constants
!USE mod_mpi

IMPLICIT NONE

CHARACTER(100)::nobj

INTEGER:: gp_N_obj
REAL(KIND=GP):: ppotV0,ppotd,ppotthick
REAL(KIND=GP):: gp_mass_part
REAL(KIND=GP):: gpdeltaE,gpr0obj
REAL(KIND=GP):: gpf0ext_amp
REAL(KIND=GP):: gpfrcpart(1:ndim)
INTEGER:: STATUS

!! Variables for use in particle dynamics calculations.

REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: gp_pos_obj,gp_temp_pos_obj
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: gp_vel_obj, gp_temp_vel_obj
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: gp_frc_obj
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: gp_f0ext_obj
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: k1_qgpp,k2_qgpp,k3_qgpp,k4_qgpp
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: k1_ugpp,k2_ugpp,k3_ugpp,k4_ugpp
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: gp_shortrange_obj

CONTAINS

SUBROUTINE create_particles

IMPLICIT NONE

PRINT*, 'NDIM', ndim

ALLOCATE(gp_pos_obj(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(gp_temp_pos_obj(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(gp_vel_obj(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(gp_temp_vel_obj(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(gp_frc_obj(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(gp_f0ext_obj(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(gp_shortrange_obj(1:ndim,1:gp_N_obj),STAT=STATUS)

ALLOCATE(k1_qgpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k2_qgpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k3_qgpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k4_qgpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k1_ugpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k2_ugpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k3_ugpp(1:ndim,1:gp_N_obj),STAT=STATUS)
ALLOCATE(k4_ugpp(1:ndim,1:gp_N_obj),STAT=STATUS)

END SUBROUTINE create_particles

SUBROUTINE destroy_particles

IMPLICIT NONE

DEALLOCATE(gp_pos_obj)
DEALLOCATE(gp_temp_pos_obj)
DEALLOCATE(gp_vel_obj)
DEALLOCATE(gp_temp_vel_obj)
DEALLOCATE(gp_frc_obj)
DEALLOCATE(gp_f0ext_obj)
DEALLOCATE(gp_shortrange_obj)
DEALLOCATE(k1_qgpp,k2_qgpp,k3_qgpp,k4_qgpp)
DEALLOCATE(k1_ugpp,k2_ugpp,k3_ugpp,k4_ugpp)

END SUBROUTINE destroy_particles

END MODULE mod_gpparticles

#endif

