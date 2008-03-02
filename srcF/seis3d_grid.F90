!******************************************************************************!
!*  This program generates grid and calculates the metric                     *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

#include "mod_macdrp.h"
!-----------------------------------------------------------------------------
program seis3d_grid
!-----------------------------------------------------------------------------

use constants_mod
use string_mod, only : string_conf
use para_mod
use math_mod
use mpi_mod
use nfseis_mod
use grid_mod
use io_mod

implicit none

!-----------------------------------------------------------------------------

real,allocatable :: gx(:),gy(:),gz(:),gxsin(:),gxcot(:)
character (len=SEIS_STRLEN) :: filenm

integer :: ierr
integer :: i,j,k,n1,n2,n3
integer :: n_i,n_j,n_k

!-----------------------------------------------------------------------------

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)

!-----------------------------------------------------------------------------

! calculate coordinate
n1=swmpi_globi(nx2,dims(1)-1)
n2=swmpi_globj(ny2,dims(2)-1)
n3=swmpi_globk(nz2,dims(3)-1)
allocate(gx(nx1:n1)); gx=0.0
allocate(gxsin(nx1:n1)); gxsin=0.0
allocate(gxcot(nx1:n1)); gxcot=0.0
allocate(gy(ny1:n2)); gy=0.0
allocate(gz(nz1:n3)); gz=0.0

do i=nx1,n1
   gx(i)=(i-ni1)*stepx+coord0(1)
   gxsin(i)=sin(gx(i))
   gxcot(i)=1.0_SP/tan(gx(i))
end do
do j=ny1,n2
   gy(j)=(j-nj1)*stepy+coord0(2)
end do
do k=nz1,n3
   gz(k)=(k-n3+LenFD)*stepz+coord0(3)
end do

print *, "output grid coordinate and locating ..."
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   !n_k=dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

!output grid
filenm=swmpi_rename_fnm(pnm_grid,fnm_metric)
call nfseis_grid_create_base(filenm,nx,ny,nz)
call nfseis_varput(filenm,'x',gx(ngx1:ngx2),(/1/),(/nx/),(/1/))
call nfseis_varput(filenm,'y',gy(ngy1:ngy2),(/1/),(/ny/),(/1/))
call nfseis_varput(filenm,'z',gz(ngz1:ngz2),(/1/),(/nz/),(/1/))
call nfseis_varput(filenm,'xsin',gxsin(ngx1:ngx2),(/1/),(/nx/),(/1/))
call nfseis_varput(filenm,'xcot',gxcot(ngx1:ngx2),(/1/),(/nx/),(/1/))

end do
end do
end do

!-----------------------------------------------------------------------------
!contains
!-----------------------------------------------------------------------------

end program seis3d_grid

