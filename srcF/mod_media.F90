!******************************************************************************!
!*  This module contains the variables for 3D media                           *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
module media_mod
!-----------------------------------------------------------------------------

use constants_mod, only : SEIS_STRLEN,SEIS_GEO,SEIS_ZERO,PI
use string_mod, only : string_conf
use math_mod
use para_mod
use mpi_mod, only :       &
    swmpi_rename_fnm
use nfseis_mod, only :    &
    nfseis_varget

implicit none
private
public :: media_fnm_init, &
          media_destroy,  &
          media_alloc,    &
          media_import

!-----------------------------------------------------------------------------

real(SP),dimension(:,:,:),allocatable,public :: &
     rho,mu,lambda,Qs
character (len=SEIS_STRLEN),public ::       &
     pnm_media,                             &
     fnm_media_conf,                        &
     fnm_media
integer :: ierr

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

subroutine media_fnm_init(fnm_conf)
  character (len=*) :: fnm_conf
  integer fid
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
    call string_conf(fid,1,'MEDIA_CONF',2,fnm_media_conf)
    call string_conf(fid,1,'MEDIA_ROOT',2,pnm_media)
    fnm_media='media.nc'
  close(fid)
end subroutine media_fnm_init

!*************************************************************************
!*                    PART-I  alloc and dealloc                          *
!*************************************************************************
subroutine media_alloc
  allocate( mu(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  mu=0.0_SP
  allocate( lambda(nx1:nx2,ny1:ny,nz1:nz2),stat=ierr);  lambda=0.0_SP
  allocate( rho(nx1:nx2,ny1:ny,nz1:nz2),stat=ierr);  rho=0.0_SP
#ifdef WITHQS
  allocate( Qs(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Qs=0.0_SP
#endif
end subroutine media_alloc
subroutine media_destroy
  if (allocated(mu)) deallocate( mu )
  if (allocated(lambda)) deallocate( lambda )
  if (allocated(rho)) deallocate( rho )
  if (allocated(Qs)) deallocate( Qs )
end subroutine media_destroy

!*************************************************************************
!*                    PART-II  media                       *
!*************************************************************************
subroutine media_import
character (len=SEIS_STRLEN) :: filenm
integer,dimension(SEIS_GEO) :: subs,subc,subt
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=swmpi_rename_fnm(pnm_media,fnm_media)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt)
  call nfseis_varget( filenm, 'rho', rho, subs,subc,subt)
#ifdef WITHQS
  call nfseis_varget( filenm, 'Qs', Qs, subs,subc,subt)
#endif
end subroutine media_import
  
end module media_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
