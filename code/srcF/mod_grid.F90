!******************************************************************************!
!*  This module contains the variables for grid geometry                      *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
module grid_mod
!-----------------------------------------------------------------------------

use constants_mod
use math_mod
use string_mod
use para_mod
use mpi_mod
use nfseis_mod

implicit none

private
public ::        &
  grid_fnm_init, &
  grid_fnm_get,  &
  grid_alloc,    &
  grid_import,   &
  grid_dealloc

!-----------------------------------------------------------------------------
real(SP),dimension(:),allocatable,public :: &
     x,y,z,xsin,xcot,                       &
     xi_x,eta_y,zeta_z
character (len=SEIS_STRLEN),public :: &
     fnm_grid_conf, pnm_grid

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

!*************************************************************************
!*                            alloc and dealloc                          *
!*************************************************************************
subroutine grid_alloc
integer :: ierr
  allocate( x(nx1:nx2),stat=ierr);  x=0.0_SP
  allocate( xsin(nx1:nx2),stat=ierr);  xsin=0.0_SP
  allocate( xcot(nx1:nx2),stat=ierr);  xcot=0.0_SP
  allocate( y(ny1:ny2),stat=ierr);  y=0.0_SP
  allocate( z(nz1:nz2),stat=ierr);  z=0.0_SP
  allocate( xi_x(nx1:nx2),stat=ierr);  xi_x=0.0_SP
  allocate( eta_y(ny1:ny2),stat=ierr);  eta_y=0.0_SP
  allocate( zeta_z(nz1:nz2),stat=ierr);  zeta_z=0.0_SP
end subroutine grid_alloc

subroutine grid_dealloc
  if (allocated(x    )) deallocate(x    )
  if (allocated(y    )) deallocate(y    )
  if (allocated(xsin )) deallocate(xsin )
  if (allocated(xcot )) deallocate(xcot )
  if (allocated(z    )) deallocate(z    )
  if (allocated(xi_x )) deallocate(xi_x )
  if (allocated(eta_y)) deallocate(eta_y)
  if (allocated(zeta_z)) deallocate(zeta_z)
end subroutine grid_dealloc

!*************************************************************************
!*                                  grid io                              *
!*************************************************************************
subroutine grid_fnm_init(fnm_conf)
character (len=*) :: fnm_conf
integer fid
fid=1001
open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'GRID_CONF',2,fnm_grid_conf)
  call string_conf(fid,1,'GRID_ROOT',2,pnm_grid)
close(fid)
end subroutine grid_fnm_init

subroutine grid_import(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=grid_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'x', x, (/1/),(/nx/),(/1/))
  call nfseis_varget( filenm, 'xsin', xsin, (/1/),(/nx/),(/1/))
  call nfseis_varget( filenm, 'xcot', xcot, (/1/),(/nx/),(/1/))
  call nfseis_varget( filenm, 'y', y, (/1/),(/ny/),(/1/))
  call nfseis_varget( filenm, 'z', z, (/1/),(/nz/),(/1/))
  call nfseis_varget( filenm, 'xi_x',   xi_x, (/1/),(/nx/),(/1/))
  call nfseis_varget( filenm, 'eta_y',  eta_y, (/1/),(/ny/),(/1/))
  call nfseis_varget( filenm, 'zeta_z', zeta_z, (/1/),(/nz/),(/1/))
end subroutine grid_import

function grid_fnm_get(n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  !filenm=trim(pnm_grid)//'/'//'coord'//'_'//set_mpi_subfix(n_i,n_j,n_k)//'.nc'
  filenm=trim(pnm_grid)//'/'//set_mpi_prefix(n_i,n_j,n_k)//'_coord.nc'
end function grid_fnm_get

end module grid_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
