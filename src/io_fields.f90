SUBROUTINE io_fields(iopoint,time)

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
!INTEGER, INTENT(IN):: jtime
INTEGER, INTENT(IN):: iopoint
REAL(KIND=GP), INTENT(IN):: time
REAL(KIND=GP):: x,y
REAL(KIND=GP):: kx,ky
REAL(KIND=GP):: rk,rk2
REAL(KIND=GP):: atemp,amin,amax,tmax,tmin
CHARACTER(LEN=100):: filename

WRITE(fnum,'(i8)') iopoint
WRITE(prc,'(i8)') (myrank+1)

!! Random number seed
CALL RANDOM_SEED(GET=randf_seed)

OPEN(UNIT=1,FILE='data/ranf_seeds/seed'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED',STATUS='UNKNOWN',IOSTAT=ios)
DO i1=1,seedsize_n
WRITE(1) randf_seed(i1)
ENDDO
CLOSE(1)

!!++++++++++++++++++++++++++++++++

#ifdef SOL_GPE
#ifdef IOEVERYPROC
OPEN(UNIT=18,FILE='data/wf/wf'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED', status='NEW', IOSTAT=ios)
WRITE(18) iopoint,time,dt
DO i2_loc=1,local_Ny_cmpx
WRITE(18) (psi(i1,i2_loc),i1=1,Nx)
ENDDO
CLOSE(18)
#endif

#ifdef IOHDF5
filename = 'data/wf/wf_'//TRIM(ADJUSTL(fnum))//'.h5'
! Separate psi into real and imaginary parts.
DO i2_loc = 1, local_Ny_cmpx
DO i1 = 1, Nx
psire(i1,i2_loc) =  REAL(psi(i1,i2_loc))
psiim(i1,i2_loc) = AIMAG(psi(i1,i2_loc))
ENDDO
ENDDO
CALL DUMP_GPE_FIELD(filename,iopoint,time)
#endif

#endif

END SUBROUTINE io_fields

!!=========================
#ifdef SOL_GPE

SUBROUTINE restart_gpe(nrestart)

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
REAL(KIND=GP):: dt_read
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

!! Read (restart from previous) random numbers
WRITE(prc,'(i8)') myrank+1
OPEN(UNIT=1,FILE='data/ranf_seeds/seed'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ios)
DO i1=1,seedsize_n
READ(1) randf_seed(i1)
ENDDO
CLOSE(1)
CALL RANDOM_SEED(PUT=randf_seed)


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

READ(11) icpoint,time,dt_read

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

READ(11) icpoint,time,dt_read

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

READ(11) icpoint,time,dt_read

DO i2_loc=1,local_Ny_cmpx
READ(11) (psi(i1,i2_loc),i1=1,Nx)
ENDDO
CLOSE(11)

ENDIF

#endif
! Closed IOEVERYPROC

#ifdef IOHDF5
filename = 'restart/wf/wf_'//TRIM(ADJUSTL(fnum))//'.h5'
CALL LOAD_GPE_FIELD(filename,icpoint,time,dt_read)
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
#ifdef DT_VAL_FROM_FILE
dt = dt_read
#endif
CALL fft(plan_psi%pname,forward)

IF (myrank .EQ. 0) THEN
PRINT*,'Restart file:',nrestart
PRINT*,'Restart no. as inside file:',icpoint
PRINT*,'Restart time:',time
PRINT*,'Time step as inside file:',dt_read
ENDIF


END SUBROUTINE restart_gpe
#endif

#ifdef IOHDF5

!!=========================
#ifdef SOL_GPE
!! Dump GPE fields
SUBROUTINE DUMP_GPE_FIELD(filename,iopoint,time)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
USE HDF5
    
IMPLICIT NONE

INTEGER, INTENT(IN):: iopoint
REAL(KIND=GP), INTENT(IN):: time

CHARACTER(LEN=100), INTENT(IN) :: filename != "wf.h5"
CHARACTER(LEN=13), PARAMETER :: dsetname1  = "psire"
CHARACTER(LEN=13), PARAMETER :: dsetname2  = "psiim"
CHARACTER(LEN=16), PARAMETER :: attrname1  = "iopoint"
CHARACTER(LEN=16), PARAMETER :: attrname2  = "time_dt"

INTEGER(HID_T):: file_id      ! File identifier
INTEGER(HID_T):: group_id     ! Group identifier
INTEGER(HID_T):: dset1_id     ! Dataset 1 identifier
INTEGER(HID_T):: dset2_id     ! Dataset 2 identifier 
INTEGER(HID_T):: filespace    ! Dataspce identifier in file
INTEGER(HID_T):: memspace     ! Dataspace identifier in memory
INTEGER(HID_T):: plist_id     ! Property list identifier

INTEGER(HID_T) :: attr1_id     ! Attribute identifier
INTEGER(HID_T) :: attr2_id     ! Attribute identifier
INTEGER(HID_T) :: aspace1_id   ! Attribute database identifier
INTEGER(HID_T) :: aspace2_id   ! Attribute database identifier
INTEGER(HID_T) :: atype1_id    ! 
INTEGER(HID_T) :: atype2_id    !

INTEGER(HSIZE_T), DIMENSION(1) :: adims1=(/1/) ! Attribute dimension
INTEGER(HSIZE_T), DIMENSION(1) :: adims2=(/2/) ! Attribute dimension
INTEGER :: arank = 1          ! Attribute rank
INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string

INTEGER, DIMENSION(1) :: attr1_data
REAL(KIND=GP), DIMENSION(1:2):: attr2_data

INTEGER(HSIZE_T), DIMENSION(2):: dimsf  != (/Nx,Ny,Nz/)
INTEGER(HSIZE_T), DIMENSION(2):: dimsfi != (/Nx,Ny,Nz/)
INTEGER(HSIZE_T), DIMENSION(2):: dimsm  != (/Nx,Ny,Nz/)
INTEGER(HSIZE_T), DIMENSION(2):: dimsmi != (/Nx,Ny,Nz/)
INTEGER:: rank = 2
INTEGER:: error

INTEGER(HSIZE_T), DIMENSION(2):: countf,  countm
INTEGER(HSIZE_T), DIMENSION(2):: offsetf, offsetm
INTEGER(HSIZE_T), DIMENSION(2):: stridef, stridem
INTEGER(HSIZE_T), DIMENSION(2):: blockf,  blockm

dimsf  = (/Nx,Ny/)
dimsfi = (/Nx,Ny/)
dimsm  = (/Nx,Ny/)
dimsmi = (/Nx,Ny/)

!WRITE(fnum,'(i8)') iopoint
!filename = 'vel_'//TRIM(ADJUSTL(fnum))//'.h5'

! Prepare psire and psiim fields.


! Set the attribute data.
attr1_data(1) = iopoint
attr2_data(1) = time
attr2_data(2) = dt

! Initialize FORTRAN predefined datatypes
CALL h5open_f(error)

! Setup file access property list for MPI-IO access.
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

! Create the file collectively.
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, &
    access_prp = plist_id)
CALL h5pclose_f(plist_id, error)

! Create scalar data space for the attribute
CALL h5screate_simple_f(arank, adims1, aspace1_id, error)
CALL h5screate_simple_f(arank, adims2, aspace2_id, error)
CALL h5acreate_f(file_id, 'Time: iopoint', H5T_NATIVE_INTEGER, aspace1_id, attr1_id, error)
CALL h5acreate_f(file_id, 'Time: time dt', H5T_NATIVE_DOUBLE, aspace2_id, attr2_id, error)

CALL h5awrite_f(attr1_id, H5T_NATIVE_INTEGER, attr1_data, adims1, error)
CALL h5awrite_f(attr2_id, H5T_NATIVE_DOUBLE,  attr2_data, adims2, error)
CALL h5aclose_f(attr1_id, error)
CALL h5aclose_f(attr2_id, error)
CALL h5sclose_f(aspace1_id, error)
CALL h5sclose_f(aspace2_id, error)

! Create the data space for the dataset.
CALL h5screate_simple_f(rank, dimsf, filespace, error)
CALL h5screate_simple_f(rank, dimsm, memspace,  error)

! Create a dataset with default properties.
CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, filespace, &
     dset1_id, error)
CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, filespace, &
     dset2_id, error)

! Each process defines dataset in memory and writes it to the
! hyperslab in the file.
! Memory hyperslab selectioon.
offsetm(1) = 0
offsetm(2) = 0
stridem(1) = 1
stridem(2) = 1
blockm(1)  = Nx
blockm(2)  = local_Ny_cmpx
countm(1)  = 1
countm(2)  = 1
     
CALL h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,offsetm,&
     countm,error,stridem,blockm)
     
! File hyperslab selection
offsetf(1) = 0
offsetf(2) = myrank*local_Ny_cmpx
stridef(1) = 1
stridef(2) = 1
blockf(1)  = Nx
blockf(2)  = local_Ny_cmpx
countf(1)  = 1
countf(2)  = 1
     
CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offsetf,&
     countf,error,stridef,blockf)

! Create property list for collective dataset write.
CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

! Write
CALL h5dwrite_f(dset1_id, H5T_NATIVE_DOUBLE, psire, dimsfi, error, &
     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

! Write
CALL h5dwrite_f(dset2_id, H5T_NATIVE_DOUBLE, psiim, dimsfi, error, &
     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

! Close the fiest dataset.
CALL h5dclose_f(dset1_id, error)
CALL h5dclose_f(dset2_id, error)

! Close the dataspace for the first dataset.
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace,error)
CALL h5pclose_f(plist_id,error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

END SUBROUTINE DUMP_GPE_FIELD

!! Load GPE fields
SUBROUTINE LOAD_GPE_FIELD(filename,icpoint,time,ddt)

USE FFTW3
USE mod_precision
USE mod_grid
USE mod_fft
USE mod_constants
USE mod_2dHD
USE mod_mpi
USE HDF5
    
IMPLICIT NONE

INTEGER, INTENT(OUT):: icpoint
REAL(KIND=GP), INTENT(OUT):: time
REAL(KIND=GP), INTENT(OUT):: ddt

CHARACTER(LEN=100), INTENT(IN) :: filename != "wf.h5"
CHARACTER(LEN=13), PARAMETER :: dsetname1  = "psire"
CHARACTER(LEN=13), PARAMETER :: dsetname2  = "psiim"
CHARACTER(LEN=16), PARAMETER :: attrname1  = "iopoint"
CHARACTER(LEN=16), PARAMETER :: attrname2  = "time_dt"

INTEGER(HID_T):: file_id      ! File identifier
INTEGER(HID_T):: group_id     ! Group identifier
INTEGER(HID_T):: dset1_id     ! Dataset 1 identifier
INTEGER(HID_T):: dset2_id     ! Dataset 2 identifier 
INTEGER(HID_T):: filespace    ! Dataspce identifier in file
INTEGER(HID_T):: memspace     ! Dataspace identifier in memory
INTEGER(HID_T):: plist_id     ! Property list identifier

INTEGER(HID_T) :: attr1_id     ! Attribute identifier
INTEGER(HID_T) :: attr2_id     ! Attribute identifier
INTEGER(HID_T) :: aspace1_id   ! Attribute database identifier
INTEGER(HID_T) :: aspace2_id   ! Attribute database identifier
INTEGER(HID_T) :: atype1_id    ! 
INTEGER(HID_T) :: atype2_id    !

INTEGER(HSIZE_T), DIMENSION(1) :: adims1=(/1/) ! Attribute dimension
INTEGER(HSIZE_T), DIMENSION(1) :: adims2=(/2/) ! Attribute dimension
INTEGER :: arank = 1          ! Attribute rank
INTEGER(SIZE_T) :: attrlen    ! Length of the attribute string

INTEGER, DIMENSION(1) :: attr1_data
REAL(KIND=GP), DIMENSION(1:2):: attr2_data

INTEGER(HSIZE_T), DIMENSION(2):: dimsf  != (/Nx,Ny,Nz/)
INTEGER(HSIZE_T), DIMENSION(2):: dimsfi != (/Nx,Ny,Nz/)
INTEGER(HSIZE_T), DIMENSION(2):: dimsm  != (/Nx,Ny,Nz/)
INTEGER(HSIZE_T), DIMENSION(2):: dimsmi != (/Nx,Ny,Nz/)
INTEGER:: rank = 2
INTEGER:: error

INTEGER(HSIZE_T), DIMENSION(2):: countf,  countm
INTEGER(HSIZE_T), DIMENSION(2):: offsetf, offsetm
INTEGER(HSIZE_T), DIMENSION(2):: stridef, stridem
INTEGER(HSIZE_T), DIMENSION(2):: blockf,  blockm

dimsf  = (/Nx,Ny/)
dimsfi = (/Nx,Ny/)
dimsm  = (/Nx,Ny/)
dimsmi = (/Nx,Ny/)

!WRITE(fnum,'(i8)') iopoint
!filename = 'vel_'//TRIM(ADJUSTL(fnum))//'.h5'

! Prepare psire and psiim fields.


! Set the attribute data.
! They are read here from the file, below.

! Initialize FORTRAN predefined datatypes
CALL h5open_f(error)

! Setup file access property list for MPI-IO access.
CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

! Create the file collectively.
CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error, &
    access_prp = plist_id)
CALL h5pclose_f(plist_id, error)

! Create scalar data space for the attribute
CALL h5screate_simple_f(arank, adims1, aspace1_id, error)
CALL h5screate_simple_f(arank, adims2, aspace2_id, error)
CALL h5aopen_f(file_id, 'Time: iopoint', attr1_id, error)
CALL h5aopen_f(file_id, 'Time: time dt', attr2_id, error)

CALL h5aread_f(attr1_id, H5T_NATIVE_INTEGER, attr1_data, adims1, error)
CALL h5aread_f(attr2_id, H5T_NATIVE_DOUBLE,  attr2_data, adims2, error)
CALL h5aclose_f(attr1_id, error)
CALL h5aclose_f(attr2_id, error)
CALL h5sclose_f(aspace1_id, error)
CALL h5sclose_f(aspace2_id, error)

! Create the data space for the dataset.
CALL h5screate_simple_f(rank, dimsf, filespace, error)
CALL h5screate_simple_f(rank, dimsm, memspace,  error)

! Open a dataset with default properties.
CALL h5dopen_f(file_id, dsetname1, dset1_id, error)
CALL h5dopen_f(file_id, dsetname2, dset2_id, error)

! Each process defines dataset in memory and writes it to the
! hyperslab in the file.
! Memory hyperslab selectioon.
offsetm(1) = 0
offsetm(2) = 0
stridem(1) = 1
stridem(2) = 1
blockm(1)  = Nx
blockm(2)  = local_Ny_cmpx
countm(1)  = 1
countm(2)  = 1
     
CALL h5sselect_hyperslab_f(memspace,H5S_SELECT_SET_F,offsetm,&
     countm,error,stridem,blockm)
     
! File hyperslab selection
offsetf(1) = 0
offsetf(2) = myrank*local_Ny_cmpx
stridef(1) = 1
stridef(2) = 1
blockf(1)  = Nx
blockf(2)  = local_Ny_cmpx
countf(1)  = 1
countf(2)  = 1
     
CALL h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,offsetf,&
     countf,error,stridef,blockf)

! Create property list for collective dataset write.
CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

! Write
CALL h5dread_f(dset1_id, H5T_NATIVE_DOUBLE, psire, dimsfi, error, &
     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

! Write
CALL h5dread_f(dset2_id, H5T_NATIVE_DOUBLE, psiim, dimsfi, error, &
     file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)

! Close the fiest dataset.
CALL h5dclose_f(dset1_id, error)
CALL h5dclose_f(dset2_id, error)

! Close the dataspace for the first dataset.
CALL h5sclose_f(filespace, error)
CALL h5sclose_f(memspace,error)
CALL h5pclose_f(plist_id,error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)

! Set the stamp/attribute information for passing back.
icpoint = attr1_data(1)
time    = attr2_data(1)
ddt     = attr2_data(2)

END SUBROUTINE LOAD_GPE_FIELD

#endif
! GPE

#endif
