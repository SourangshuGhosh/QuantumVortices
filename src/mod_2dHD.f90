MODULE mod_2dHD

USE mod_precision
USE FFTW3
USE mod_grid

IMPLICIT NONE

INTEGER:: ioerr
CHARACTER(256):: iomesg

!--------
! Random number-seed
INTEGER, ALLOCATABLE :: randf_seed(:)
INTEGER :: seedsize_n
INTEGER :: minitseed

!!------------GPE------------------------
#ifdef SOL_GPE
REAL(KIND=GP):: Nparticle
REAL(KIND=GP):: gstrength
REAL(KIND=GP):: csound
REAL(KIND=GP):: xi
REAL(KIND=GP):: alpha
REAL(KIND=GP):: omega_chem
REAL(KIND=GP):: rho_avg
REAL(KIND=GP):: noisesigma
REAL(KIND=GP):: oneybeta
REAL(KIND=GP):: nuN
REAL(KIND=GP):: chempot0
REAL(KIND=GP):: disp_factor
#ifdef GPE_COUNTERFLOW
REAL(KIND=GP):: vnx,vny
#endif

#ifdef GPE_SELFTRUNC
REAL(KIND=GP):: ktrunc
#endif

REAL(KIND=GP):: nk0
INTEGER, ALLOCATABLE, DIMENSION(:):: kwn1_c2c,kwn2_c2c
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: spec_gp

!! Forcing and dissipation. Wave turbulence.
#ifdef GPE_FRCDISP
REAL(KIND=GP):: famp_gp, kf1_gp, kf2_gp
INTEGER:: nhypvis1_gp, nhypvis2_gp
REAL(KIND=GP):: vishypo_gp, vishype_gp
REAL(KIND=GP):: friclowk_gp
#endif

#ifdef GPOTOC
REAL(KIND=GP):: k_otoc_gp,eps_otoc_gp
#endif

#endif
!!--------------------------------------


CONTAINS

SUBROUTINE create_global_array

IMPLICIT NONE

CALL RANDOM_SEED(SIZE = seedsize_n)
ALLOCATE(randf_seed(1:seedsize_n))

!!---------------------------------------
#ifdef SOL_GPE
ALLOCATE(kwn1_c2c(1:Nx))
ALLOCATE(kwn2_c2c(1:Ny))
ALLOCATE(spec_gp(0:nshell,1:6))
#endif

END SUBROUTINE create_global_array

SUBROUTINE dealloc_global_array

IMPLICIT NONE


DEALLOCATE(randf_seed)

#ifdef SOL_GPE
DEALLOCATE(kwn1_c2c,kwn2_c2c)
DEALLOCATE(spec_gp)
#endif
!!---------------------------------------

END SUBROUTINE dealloc_global_array

END MODULE mod_2dHD
