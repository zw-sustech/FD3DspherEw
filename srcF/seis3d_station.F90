!******************************************************************************!
!*  This program set distribute station among mpi threads                     *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
program seis3d_station
!-----------------------------------------------------------------------------

#include "mod_macdrp.h"

use constants_mod
use string_mod, only : string_conf
use para_mod
use math_mod
use mpi_mod
use nfseis_mod
use grid_mod
use io_mod

implicit none

integer :: n_i,n_j,n_k,n1,n2,n3
real,allocatable :: gx(:),gy(:),gz(:)

!-----------------------------------------------------------------------------

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)
call grid_alloc
call io_init(fnm_conf)
call io_pt_read(fnm_conf)

call grid_alloc

n1=swmpi_globi(nx2,dims(1)-1)
n2=swmpi_globj(ny2,dims(2)-1)
n3=swmpi_globk(nz2,dims(3)-1)
allocate(gx(nx1:n1)); gx=0.0
allocate(gy(ny1:n2)); gy=0.0
allocate(gz(nz1:n3)); gz=0.0

do n_i=0,dims(1)-1
   n_j=0; n_k=0
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_import
   gx(ngx1:ngx2)=x
end do

do n_j=0,dims(2)-1
   n_i=0; n_k=0
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_import
   gy(ngy1:ngy2)=y
end do

do n_k=0,dims(3)-1
   n_i=0; n_j=0
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_import
   gz(ngz1:ngz2)=z
end do

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1

  write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
  call swmpi_change_fnm(n_i,n_j,n_k)
  call swmpi_set_gindx(n_i,n_j,n_k)

  call receiver_locate(n_i,n_j,n_k)

end do
end do
end do

call grid_dealloc

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine receiver_locate(n_i,n_j,n_k)
integer,intent(in) :: n_i,n_j,n_k
integer,dimension(SEIS_GEO) :: p1
character (len=SEIS_STRLEN) :: filenm
integer npt,n,i,j,k,gi,gj,gk
real(SP) x0,y0,z0
npt=0
filenm=get_fnm_seismoinfo(pnm_seismoinfo,n_i,n_j,n_k)
call nfseis_seismoinfo_create(filenm,pt_tinv)
do n=1,num_pt
   x0=pt_xyz(1,n);y0=pt_xyz(2,n);z0=pt_xyz(3,n)
   p1(1:1)=minloc(abs(gx-x0)); p1(2:2)=minloc(abs(gy-y0))
   gi=loct_i(p1(1));gj=loct_j(p1(2))
if (      gi>=ngi1 .and. gi<=ngi2 &
    .and. gj>=ngj1 .and. gj<=ngj2 ) then
   i=swmpi_locli(gi,n_i)
   j=swmpi_loclj(gj,n_j)
   if (z0<=topo_hyper_height) then
      p1(3:3)=minloc(abs(gz-z0))
      gk=loct_k(p1(3))
      k=swmpi_loclk(gk,n_k)
   elseif (n_k==dims(3)-1) then
      gk=ngk2; k=swmpi_loclk(gk,n_k)
      z0=gz(gk)
   end if
   if (gk>=ngk1 .and. gk<=ngk2) then
      npt=npt+1
      gi=swmpi_globi(i,n_i)
      gj=swmpi_globj(j,n_j)
      gk=swmpi_globk(k,n_k)
      call nfseis_varput(filenm,'indx',(/out_i(i),out_j(j),out_k(k)/),     &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'gindx',(/out_i(gi),out_j(gj),out_k(gk)/), &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'coord',(/x0,y0,z0/),                      &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'grid',(/x(i),y(j),z(k)/),                 &
           (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'id',pt_id(:,n),                           &
           (/1,npt/),(/2,1/),(/1,1/))
   end if
end if
end do
end subroutine receiver_locate

end program seis3d_station

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
