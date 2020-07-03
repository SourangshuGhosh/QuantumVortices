MODULE mod_mpi

IMPLICIT NONE
SAVE

INCLUDE "mpif.h"

INTEGER:: istatus(MPI_STATUS_SIZE)
INTEGER:: ierr,nprocs,myrank

END MODULE mod_mpi
