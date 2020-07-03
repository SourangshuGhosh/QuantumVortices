MODULE mod_fft

USE FFTW3
USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD
USE mod_mpi

IMPLICIT NONE

INTEGER:: local_Ny,local_i2_offset
INTEGER:: local_Nx,local_i1_offset
!! FFTW definitions make use of C_INTPTR/C_INT.
!! We copy Nx,Ny,Nz into suitable data types.
INTEGER(C_INTPTR_T):: RNx
INTEGER(C_INTPTR_T):: RNy

!INTEGER(C_INTPTR_T):: RNxh,RNyh,RNzh
!INTEGER(C_INTPTR_T):: RNxhp,RNyhp,RNzhp
!INTEGER(C_INTPTR_T):: RNxpp,RNypp,RNzpp

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

!!++++++++++++++++++++++++++++++++++++++++++

!--------------GPE-----------------
#ifdef SOL_GPE
TYPE(plan):: plan_psi
TYPE(plan):: plan_gpnlin
TYPE(plan):: plan_cmpxtmp
TYPE(plan):: plan_gptmp1
TYPE(plan):: plan_gptmp2
TYPE(plan):: plan_gptmp3
#ifdef GPPART
TYPE(plan):: plan_gptransvobj
#endif
#endif
!!----------------------------------


!!++++++++++++++++++++++++++++++++++++++++++
!!--------------BEGIN---GPE-----------------
#ifdef SOL_GPE
TYPE(C_PTR):: data_psi,data_kpsi
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: psi(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: kpsi(:,:)

TYPE(C_PTR):: data_gpnlin
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_nlin(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_nlin(:,:)

TYPE(C_PTR):: data_cmpxtmp
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: phys_cmpxtmp(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: four_cmpxtmp(:,:)
#endif

!!++++++++++++++++++++++++++++++++++++++++++
!!--------------BEGIN---GPE-----------------
#ifdef SOL_GPE
TYPE(C_PTR):: data_gpstr1
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_str1(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_str1(:,:)
REAL(C_DOUBLE), POINTER:: psire(:,:)

TYPE(C_PTR):: data_gpstr2
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_str2(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_str2(:,:)
REAL(C_DOUBLE), POINTER:: psiim(:,:)

TYPE(C_PTR):: data_gptmp1
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_tmp1(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_tmp1(:,:)

TYPE(C_PTR):: data_gptmp2
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_tmp2(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_tmp2(:,:)

TYPE(C_PTR):: data_gptmp3
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_tmp3(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_tmp3(:,:)

#ifdef GPPART
TYPE(C_PTR):: data_vobj,data_trans_vobj
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_vobj(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_vobj(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_phys_trans_vobj(:,:)
COMPLEX(C_DOUBLE_COMPLEX), POINTER:: cgp_four_trans_vobj(:,:)
#endif
#endif

CONTAINS

SUBROUTINE create_state_variables

IMPLICIT NONE

dimnd=2

IF (myrank==0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*,'Start Allocation of vel, etc.'
PRINT*, '!!---------------'
ENDIF

!! Memory allocation and pointer association 
!! for state variables.

!!---------------------------------
! GPE ALLOCATION
!!---------------------------------
#ifdef SOL_GPE
!!-- Complex wave function.
data_psi = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_psi,psi,[Nx,local_Ny_cmpx])

data_kpsi = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_kpsi,kpsi,[Ny,local_Nx_cmpx])

data_gpnlin = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_gpnlin,cgp_phys_nlin,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_gpnlin,cgp_four_nlin,[Ny,local_Nx_cmpx])

data_cmpxtmp = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_cmpxtmp,phys_cmpxtmp,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_cmpxtmp,four_cmpxtmp,[Ny,local_Nx_cmpx])

data_gpstr1 = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_gpstr1,cgp_phys_str1,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_gpstr1,cgp_four_str1,[Ny,local_Nx_cmpx])
!! psire: only for dumping real part of psi.
CALL c_f_pointer(data_gpstr1,psire,[Ny,2*local_Nx_cmpx])

data_gpstr2 = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_gpstr2,cgp_phys_str2,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_gpstr2,cgp_four_str2,[Ny,local_Nx_cmpx])
!! psiim: only for dumping imaginary part of psi.
CALL c_f_pointer(data_gpstr2,psiim,[Ny,2*local_Nx_cmpx])

data_gptmp1 = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_gptmp1,cgp_phys_tmp1,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_gptmp1,cgp_four_tmp1,[Ny,local_Nx_cmpx])

data_gptmp2 = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_gptmp2,cgp_phys_tmp2,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_gptmp2,cgp_four_tmp2,[Ny,local_Nx_cmpx])

data_gptmp3 = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_gptmp3,cgp_phys_tmp3,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_gptmp3,cgp_four_tmp3,[Ny,local_Nx_cmpx])

#ifdef GPPART
data_vobj = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_vobj,cgp_phys_vobj,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_vobj,cgp_four_vobj,[Ny,local_Nx_cmpx])

data_trans_vobj = fftw_alloc_complex(alloc_local_cmpx)
CALL c_f_pointer(data_trans_vobj,cgp_phys_trans_vobj,[Nx,local_Ny_cmpx])
CALL c_f_pointer(data_trans_vobj,cgp_four_trans_vobj,[Ny,local_Nx_cmpx])
#endif

#endif

IF (myrank==0) THEN 
PRINT*, '!!-----------------------------------------------!!'
PRINT*,'Vel etc. allocated'
PRINT*, '!!-----------------------------------------------!!'
ENDIF

END SUBROUTINE create_state_variables


SUBROUTINE free_state_variables

IMPLICIT NONE

!! Release the allocated memory.

!!---------------------------------
#ifdef SOL_GPE
CALL fftw_free(data_psi)
CALL fftw_free(data_kpsi)
CALL fftw_free(data_gpnlin)
CALL fftw_free(data_cmpxtmp)

CALL fftw_free(data_gpstr1)
CALL fftw_free(data_gpstr2)

CALL fftw_free(data_gptmp1)
CALL fftw_free(data_gptmp2)
CALL fftw_free(data_gptmp3)

#ifdef GPPART
CALL fftw_free(data_vobj)
CALL fftw_free(data_trans_vobj)
#endif

#endif

END SUBROUTINE free_state_variables

SUBROUTINE mpi_localsize_estimate

IMPLICIT NONE

IF (myrank==0) THEN
PRINT*, '!!-----------------------------------------------!!'
PRINT*, '!! Estimation of local array sizes using FFTW.'
PRINT*, 'Use: fftw_mpi_local_size_transposed'
ENDIF

!======================================================================
!       Definitions for the fftw operations
!======================================================================

!! Obtain the size of the local arrays on each process, using FFTW
!! slab decomposition.
!!--------------------------------------------
RNx=Nx
RNy=Ny

! --- Get the local size for 3d arrays.
lclsize2%rank = 2
lclsize2%nrank(1:2) = (/RNy,RNx/)

#ifdef SOL_GPE
!alloc_local_cmpx = fftw_mpi_local_size_3d_transposed(Nz,Ny,Nx,&
!        MPI_COMM_WORLD,local_Nz_cmpx,local_i3_offset_cmpx,&
!                       local_Ny_cmpx,local_i2_offset_cmpx)
alloc_local_cmpx = fftw_mpi_local_size_transposed(lclsize2%rank,lclsize2%nrank,&
MPI_COMM_WORLD,Rlocal_Ny,Rlocal_i2_offset,Rlocal_Nx,Rlocal_i1_offset)
local_Ny_cmpx=Rlocal_Ny
local_Nx_cmpx=Rlocal_Nx
local_i2_offset_cmpx=Rlocal_i2_offset
local_i1_offset_cmpx=Rlocal_i1_offset

IF (myrank .EQ. 0) THEN
!!----------------------
PRINT*, 'For complex to complex transforms'
!!----------------------
PRINT*, '!! Physical space'
PRINT*, 'Size when slabbed along y-direction (alloc_cmpx):', alloc_local_cmpx
PRINT*, 'Local_Nz:', local_Ny_cmpx
PRINT*, 'local_i3_offset:', local_i2_offset_cmpx
!!----------------------
PRINT*, '!! Fourier space'
PRINT*, 'In Fourier space the fields are left in a transposed state compared to Phys. spc.'
PRINT*, 'This means that the y-dimension is contained in the third index of the array,'
PRINT*, 'along which the array is swaped.'
PRINT*, 'Local_Ny:', local_Nx_cmpx
PRINT*, 'local_i2_offset', local_i1_offset_cmpx
!!-----------------------------------------------!!
ENDIF
#endif

IF (myrank==0) PRINT*,'Check code status after estimating localsize'

END SUBROUTINE mpi_localsize_estimate

SUBROUTINE initialize_fft

IMPLICIT NONE

!! No. of grid points copied to variables  with suitable datatypes for fftw.
RNx=Nx
RNy=Ny

!!-------------BEGIN---GPE--BLOCK----------------
#ifdef SOL_GPE
!!----
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

!!----- GP Nonlinear term.
plan_gpnlin%pname = 'gpnlin'
plan_gpnlin%pflag = 1
plan_gpnlin%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_gpnlin%nrank(1:2)=(/RNy,RNx/)
IF (myrank .EQ. 0) PRINT*, 'Plan gpnlin'
!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_gpnlin%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_phys_nlin,cgp_four_nlin,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_gpnlin%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_four_nlin,cgp_phys_nlin,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)

!!-----
plan_cmpxtmp%pname = 'cmpxtmp'
plan_cmpxtmp%pflag = 1
plan_cmpxtmp%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_cmpxtmp%nrank(1:2)=(/RNy,RNx/)
IF (myrank .EQ. 0) PRINT*, 'Plan cmpxtmp'
!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_cmpxtmp%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       phys_cmpxtmp,four_cmpxtmp,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_cmpxtmp%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       four_cmpxtmp,phys_cmpxtmp,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)

!!-----
plan_gptmp1%pname = 'gptmp1'
plan_gptmp1%pflag = 1
plan_gptmp1%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_gptmp1%nrank(1:2)=(/RNy,RNx/)
IF (myrank .EQ. 0) PRINT*, 'Plan gptmp1'
!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_gptmp1%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_phys_tmp1,cgp_four_tmp1,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_gptmp1%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_four_tmp1,cgp_phys_tmp1,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)

!!-----
plan_gptmp2%pname = 'gptmp2'
plan_gptmp2%pflag = 1
plan_gptmp2%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_gptmp2%nrank(1:2)=(/RNy,RNx/)
IF (myrank .EQ. 0) PRINT*, 'Plan gptmp2'
!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_gptmp2%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_phys_tmp2,cgp_four_tmp2,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_gptmp2%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_four_tmp2,cgp_phys_tmp2,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)

!!-----
plan_gptmp3%pname = 'gptmp3'
plan_gptmp3%pflag = 1
plan_gptmp3%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_gptmp3%nrank(1:2)=(/RNy,RNx/)
IF (myrank .EQ. 0) PRINT*, 'Plan gptmp3'
!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_gptmp3%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_phys_tmp3,cgp_four_tmp3,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_gptmp3%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_four_tmp3,cgp_phys_tmp3,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)

#ifdef GPPART
!!-----
plan_gptransvobj%pname = 'gptransvobj'
plan_gptransvobj%pflag = 1
plan_gptransvobj%rank  = 2
! Logical array size of the fft arrays, in reverse order
plan_gptransvobj%nrank(1:2)=(/RNy,RNx/)
IF (myrank .EQ. 0) PRINT*, 'Plan gptransvobj'
!Nz,Ny,Nx plan_psi%rank,plan_psi%nrank
!!*** Note for some reason, plan does not accept the rank type 
plan_gptransvobj%pfor = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_phys_trans_vobj,cgp_four_trans_vobj,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_OUT)
!
plan_gptransvobj%pinv = fftw_mpi_plan_dft_2d(RNy,RNx,&
       cgp_four_trans_vobj,cgp_phys_trans_vobj,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_MEASURE+FFTW_MPI_TRANSPOSED_IN)
#endif

#endif

IF (myrank .EQ. 0) THEN 
PRINT*, '!!-----------------------------------------------!!'
PRINT*, 'FFTW plan for force created.'
PRINT*, '!!-----------------------------------------------!!'
ENDIF

!fftwscale = one/real((Nx*Ny*Nz),kind=GP)
fftwscale = one/(REAL(Nx,KIND=GP)*REAL(Ny,KIND=GP))

IF (myrank .EQ. 0) THEN
PRINT*, 'FFTW transforms are unnormalized.'
PRINT*, 'Forward transfrom followed by an inverse transform'
PRINT*, 'will leave the array multiplied by Nx*Ny.'
PRINT*, 'FFTW normalization factor (fftwscale = one/real((Nx*Ny),kind=GP)):', fftwscale
ENDIF

IF (myrank .EQ. 0) THEN
PRINT*, 'Initialization of FFT related stuff done.'
PRINT*, '!!-----------------------------------------------!!'
ENDIF

END SUBROUTINE initialize_fft

!! Subroutine to free memory associated with fftw pointers
!! and destroy fftw plans.

SUBROUTINE clear_fft

IMPLICIT NONE

!!-----------BEGIN--GPE---BLOCK---------
#ifdef SOL_GPE
CALL fftw_destroy_plan(plan_psi%pfor)
CALL fftw_destroy_plan(plan_psi%pinv)

CALL fftw_destroy_plan(plan_gpnlin%pfor)
CALL fftw_destroy_plan(plan_gpnlin%pinv)

CALL fftw_destroy_plan(plan_cmpxtmp%pfor)
CALL fftw_destroy_plan(plan_cmpxtmp%pinv)

CALL fftw_destroy_plan(plan_gptmp1%pfor)
CALL fftw_destroy_plan(plan_gptmp1%pinv)

CALL fftw_destroy_plan(plan_gptmp2%pfor)
CALL fftw_destroy_plan(plan_gptmp2%pinv)

CALL fftw_destroy_plan(plan_gptmp3%pfor)
CALL fftw_destroy_plan(plan_gptmp3%pinv)

#ifdef GPPART
CALL fftw_destroy_plan(plan_gptransvobj%pfor)
CALL fftw_destroy_plan(plan_gptransvobj%pinv)
#endif

#endif
!!-----------END---GPE---BLOCK----------

END SUBROUTINE clear_fft

!!---------------------------------------------
! FFT Subroutine for: Physical to Fourier space
!!---------------------------------------------

SUBROUTINE fft(plan_name,direction)

IMPLICIT NONE

CHARACTER (100), INTENT(IN):: plan_name
INTEGER, INTENT(IN):: direction
!INTEGER, INTENT(IN):: trn_list
INTEGER:: plan_flag

!!----------BEGIN---GPE---BLOCK---------
#ifdef SOL_GPE
IF (plan_name .EQ. plan_psi%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_psi%pfor,psi,kpsi)
ELSEIF (plan_name .EQ. plan_psi%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_psi%pinv,kpsi,psi)
psi = fftwscale*psi
ENDIF

!!--- gpnlin
IF (plan_name .EQ. plan_gpnlin%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_gpnlin%pfor,cgp_phys_nlin,cgp_four_nlin)
ELSEIF (plan_name .EQ. plan_gpnlin%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_gpnlin%pinv,cgp_four_nlin,cgp_phys_nlin)
cgp_phys_nlin = fftwscale*cgp_phys_nlin
ENDIF

!!--- cmpxtmp
IF (plan_name .EQ. plan_cmpxtmp%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_cmpxtmp%pfor,phys_cmpxtmp,four_cmpxtmp)
ELSEIF (plan_name .EQ. plan_cmpxtmp%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_cmpxtmp%pinv,four_cmpxtmp,phys_cmpxtmp)
phys_cmpxtmp = fftwscale*phys_cmpxtmp
ENDIF

!!--- gptmp1
IF (plan_name .EQ. plan_gptmp1%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_gptmp1%pfor,cgp_phys_tmp1,cgp_four_tmp1)
ELSEIF (plan_name .EQ. plan_gptmp1%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_gptmp1%pinv,cgp_four_tmp1,cgp_phys_tmp1)
cgp_phys_tmp1 = fftwscale*cgp_phys_tmp1
ENDIF

!!--- gptmp
IF (plan_name .EQ. plan_gptmp2%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_gptmp2%pfor,cgp_phys_tmp2,cgp_four_tmp2)
ELSEIF (plan_name .EQ. plan_gptmp2%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_gptmp2%pinv,cgp_four_tmp2,cgp_phys_tmp2)
cgp_phys_tmp2 = fftwscale*cgp_phys_tmp2
ENDIF

!!--- gptmp
IF (plan_name .EQ. plan_gptmp3%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_gptmp3%pfor,cgp_phys_tmp3,cgp_four_tmp3)
ELSEIF (plan_name .EQ. plan_gptmp3%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_gptmp3%pinv,cgp_four_tmp3,cgp_phys_tmp3)
cgp_phys_tmp3 = fftwscale*cgp_phys_tmp3
ENDIF

!!--- gptransvobj
#ifdef GPPART
IF (plan_name .EQ. plan_gptransvobj%pname .AND. direction .EQ. forward) THEN
CALL fftw_mpi_execute_dft(plan_gptransvobj%pfor,cgp_phys_trans_vobj,cgp_four_trans_vobj)
ELSEIF (plan_name .EQ. plan_gptransvobj%pname .AND. direction .EQ. inverse) THEN
CALL fftw_mpi_execute_dft(plan_gptransvobj%pinv,cgp_four_trans_vobj,cgp_phys_trans_vobj)
cgp_phys_trans_vobj = fftwscale*cgp_phys_trans_vobj
ENDIF
#endif

#endif
!!----------END---GPE---BLOCK-----------

END SUBROUTINE fft

END MODULE mod_fft
