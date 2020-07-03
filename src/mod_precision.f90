module mod_precision

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!       Purpose:
!       Determines the data-type to be used in the code.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

implicit none

!Double precision
integer, parameter:: GP = KIND(0.0D0)
!Single precision
!integer, parameter:: GP = KIND(0.0)

end module mod_precision
