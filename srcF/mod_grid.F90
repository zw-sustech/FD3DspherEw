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
use string_mod, only : &
    string_conf
use para_mod
use mpi_mod, only :    &
    swmpi_rename_fnm
use nfseis_mod, only : &
    nfseis_varget

implicit none

private
public :: grid_fnm_init, &
          grid_alloc,    &
          grid_import,   &
          grid_dealloc

!-----------------------------------------------------------------------------
real(SP),dimension(:),allocatable,public :: &
     x,y,z,xsin,xcot
character (len=SEIS_STRLEN),public :: &
     fnm_grid_conf,                   &
     pnm_grid,                        &
     fnm_metric

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
subroutine grid_fnm_init(fnm_conf)
character (len=*) :: fnm_conf
integer fid
fid=1001
open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'fnm_grid_conf',2,fnm_grid_conf)
  call string_conf(fid,1,'pnm_grid',2,pnm_grid)
  call string_conf(fid,1,'fnm_metric',2,fnm_metric)
close(fid)
end subroutine grid_fnm_init

!*************************************************************************
!*                    PART-I  alloc and dealloc                          *
!*************************************************************************
subroutine grid_alloc
integer :: ierr
  allocate( x(nx1:nx2),stat=ierr);  x=0.0_SP
  allocate( xsin(nx1:nx2),stat=ierr);  xsin=0.0_SP
  allocate( xcot(nx1:nx2),stat=ierr);  xcot=0.0_SP
  allocate( y(ny1:ny2),stat=ierr);  y=0.0_SP
  allocate( z(nz1:nz2),stat=ierr);  z=0.0_SP
end subroutine grid_alloc

subroutine grid_dealloc
  if (allocated(x    )) deallocate(x    )
  if (allocated(y    )) deallocate(y    )
  if (allocated(xsin )) deallocate(xsin )
  if (allocated(xcot )) deallocate(xcot )
  if (allocated(z    )) deallocate(z    )
end subroutine grid_dealloc

!*************************************************************************
!*                    PART-II  grid import                        *
!*************************************************************************

subroutine grid_import
character (len=SEIS_STRLEN) :: filenm
filenm=swmpi_rename_fnm(pnm_grid,fnm_metric)
   call nfseis_varget( filenm, 'x', x, (/1/),(/nx/),(/1/))
   call nfseis_varget( filenm, 'xsin', xsin, (/1/),(/nx/),(/1/))
   call nfseis_varget( filenm, 'xcot', xcot, (/1/),(/nx/),(/1/))
   call nfseis_varget( filenm, 'y', y, (/1/),(/ny/),(/1/))
   call nfseis_varget( filenm, 'z', z, (/1/),(/nz/),(/1/))
end subroutine grid_import

end module grid_mod
