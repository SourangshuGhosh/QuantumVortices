SUBROUTINE solve

USE mod_precision
USE mod_constants
USE mod_grid
USE mod_2dHD
USE mod_fft
USE mtmod
#ifdef GPPART
USE mod_gpparticles
#endif

IMPLICIT NONE

!INTEGER:: itime,time,i1,iseed,ickp
INTEGER:: itime,jtime,i1,iseed,ickp
REAL(KIND=GP):: time
integer:: fldpoint,spcpoint,stspoint


!! Preservation of seed for random numbers
!! on different cores for restart is dropped
!! for the time being. We are not using
!! random numbers in forcing. We should revisit this.

!iseed = (371+myrank)*1237
!CALL sgrnd(iseed)
!IF (nrst .EQ. 0) THEN
!DO i1=1,seedsize_n
!randf_seed(i1)=12345+(i1+myrank)*1237
!ENDDO
!ENDIF

minitseed = 123456789


IF (nrst .EQ. 0) THEN
itime=0
jtime=0
time = zero
fldpoint=0
spcpoint=0
stspoint=0

!! Set the random number seeds
!CALL RANDOM_SEED(GET=randf_seed)
WRITE(fnum,'(i8)') fldpoint
WRITE(prc,'(i8)') myrank+1
OPEN(UNIT=1,FILE='data/ranf_seeds/seed'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
FORM='UNFORMATTED',STATUS='NEW',IOSTAT=ios)
DO i1=1,seedsize_n
randf_seed(i1) = (myrank+1)*minitseed*i1
WRITE(1) randf_seed(i1)
ENDDO
CLOSE(1)
CALL RANDOM_SEED(PUT=randf_seed)

!CALL init_random_seed(fldpoint,myrank)
!WRITE(fnum,'(i8)') fldpoint
!WRITE(prc,'(i8)') myrank+1
!OPEN(UNIT=1,FILE='data/ranf_seeds/seed_ini_'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
!FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ios)
!DO i1=1,seedsize_n
!READ(1) randf_seed(i1)
!ENDDO
!CLOSE(1)
!CALL RANDOM_SEED(PUT=randf_seed)

ELSE
itime=0
time = timestart
fldpoint = nrst
spcpoint = sprst
stspoint = strst

!! Read (restart from previous) random numbers
!WRITE(fnum,'(i8)') fldpoint
!WRITE(prc,'(i8)') myrank+1
!OPEN(UNIT=1,FILE='data/ranf_seeds/seed'//TRIM(ADJUSTL(fnum))//'p'//TRIM(ADJUSTL(prc))//'.dat',&
!FORM='UNFORMATTED',STATUS='OLD',IOSTAT=ios)
!DO i1=1,seedsize_n
!READ(1) randf_seed(i1)
!ENDDO
!CLOSE(1)
!CALL RANDOM_SEED(PUT=randf_seed)

ENDIF

#ifdef FIXDT
IF (MOD(itime,fsteps)==0 .AND. nrst .EQ. 0) CALL io_fields(fldpoint,time)
#ifdef SOL_GPE
IF (MOD(itime,gsteps)==0 .AND. nrst .EQ. 0) CALL diagnostics_gpe(spcpoint,itime,time)
#ifdef GPE_SGLE
IF (MOD(itime,fsteps)==0 .AND. nrst .EQ. 0) CALL write_chempot_for_restart(time,omega_chem)
#endif
#ifdef GPSTS
IF (MOD(itime,ststeps)==0 .AND. nrst .EQ.0) CALL sts_gpe(stspoint,time)
#endif
#endif
#ifdef GPPART
IF (MOD(itime,fsteps)==0 .AND. nrst .EQ. 0) CALL gppart_checkpoint(fldpoint,itime,time)
#endif
#endif

DO itime=1,maxiter

time = time + dt

#ifdef SOL_GPE
#ifdef GPE_ARGLE
CALL ARGLE(time)
#ifdef GPE_SGLE
IF (MOD(itime,fsteps)==0) CALL write_chempot_for_restart(time,omega_chem)
#endif
#else
!if (myrank==0) print*,'Hello GP'
CALL gpe_rk4
!CALL ARGLE
!CALL implicit_euler_sargle
#endif
#endif

#ifdef FIXDT
IF (MOD(itime,fsteps)==0) THEN
fldpoint = fldpoint + 1
CALL io_fields(fldpoint,time)
#ifdef GPPART
CALL gppart_checkpoint(fldpoint,itime,time)
#endif
ENDIF
IF (MOD(itime,gsteps)==0) THEN
IF (MOD(itime,ssteps)==0) spcpoint = spcpoint + 1
#ifdef SOL_GPE
CALL diagnostics_gpe(spcpoint,itime,time)
#endif
ENDIF
#ifdef GPSTS
IF (MOD(itime,ststeps)==0) THEN
stspoint = stspoint + 1
CALL sts_gpe(stspoint,time)
ENDIF
#endif
#endif


ENDDO


END SUBROUTINE solve
