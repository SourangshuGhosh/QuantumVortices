module FFTW3
use, intrinsic :: iso_c_binding
!!----------Local/personal desktops----
!include 'include/fftw3-mpi.f03'
include 'include/fftw3.f03'
!!----------Bluegene-IDRIS-Turing-------
!include '/bglocal/cn/pub/FFTW/3.3.3/include/fftw3-mpi.f03'
!include '/linkhome/rech/xhr/rxhr001/fftw/fftw-3.3.4/include/fftw3-mpi.f03'
!!---------CICADA Unice-----------------
!include '/softs/fftw/fftw-3.3.3/fftw-3.3.3_intel-12.1/include/fftw3-mpi.f03'
!!---------LICALLO OCA------------------
!include '/trinity/shared/OCA/softs/fftw-3.3.7-intel17//include/fftw3-mpi.f03'
!!---------OCCIGEN CINES----------------
!include '/opt/software/occigen/libraries/fftw/3.3.6-pl2/intel/17.0/intelmpi/2017.0.098/include/fftw3-mpi.f03'
end module FFTW3
