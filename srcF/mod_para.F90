!******************************************************************************!
!*  This module declares parameters variables                                 *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

#include "mod_macdrp.h"
!-----------------------------------------------------------------------------
module para_mod
!-----------------------------------------------------------------------------

use constants_mod
use string_mod, only : string_conf
implicit none
!private

public :: para_init,get_conf_name, &
          inn_i,inn_j,inn_k ,      &
          out_i,out_j,out_k ,      &
          loct_i,loct_j,loct_k

!-----------------------------------------------------------------------------
integer,public ::                     &
    ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &! without additional stencil points
    nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &! include boundary stencil points
    ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &! global index without ghost points
    ngx1,ngx2,ngy1,ngy2,ngz1,ngz2 ! same as above with ghost
integer,dimension(SEIS_GEO*2),public :: &
    point_in_this
real,dimension(SEIS_GEO),public ::      &
    coord0
integer,public :: nt
real(SP),public :: stept,stepx,stepy,stepz
real(SP),public :: cur_time, cur_nt
character (len=SEIS_STRLEN),public :: fnm_conf

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
subroutine get_conf_name(fnm_input)
character (*),intent(out) :: fnm_input

#ifdef USE_GETARG
integer numarg,i,iargc
   numarg = iargc( )
   do i=1,numarg
      call getarg(i,fnm_input)
      !print *, "the ",i," input arg is ", trim(fnm_input)
      print *, "the conf file name is "//trim(fnm_input)
   enddo
   if ( numarg/=1 ) then
      print *, "numarg should be equal to 1, but now it is ", numarg
      stop 1
   endif
#else
   fnm_input="SeisFD3D.conf"
#endif
end subroutine get_conf_name

!---------------------------------------------------------------------------

subroutine para_init(fnm_conf)
character (len=*) :: fnm_conf
integer fid,n
fid=1001
open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'ni',2,ni)
  call string_conf(fid,1,'nj',2,nj)
  call string_conf(fid,1,'nk',2,nk)
  call string_conf(fid,1,'nt',2,nt)
  nx1=1; ni1=nx1+LenFD; ni2=ni1+ni-1; nx2=ni2+LenFD
  ny1=1; nj1=ny1+LenFD; nj2=nj1+nj-1; ny2=nj2+LenFD
  nz1=1; nk1=nz1+LenFD; nk2=nk1+nk-1; nz2=nk2+LenFD
  nx=nx2-nx1+1; ny=ny2-ny1+1; nz=nz2-nz1+1

  call string_conf(fid,1,'stept',2,stept)
  call string_conf(fid,1,'dtheta',2,stepx)
  call string_conf(fid,1,'dphi',2,stepy)
  call string_conf(fid,1,'dr',2,stepz)
  stepx=stepx*PI/180.0_SP
  stepy=stepy*PI/180.0_SP
do n=1,SEIS_GEO
  call string_conf(fid,1,'X1Y1Ztop',n+1,coord0(n))
end do
  coord0(1:2)=coord0(1:2)*PI/180.0_SP
close(fid)
  cur_time=0.0
  cur_nt  = 0
end subroutine para_init

subroutine set_cur_time(ntime,incr)
integer :: ntime
real(SP) :: incr
cur_nt = ntime
cur_time = (ntime+incr)*stept
end subroutine set_cur_time

function inn_i(oi) result(ii)
  integer,intent(in) :: oi
  integer :: ii
  ii=oi+ni1-1
end function inn_i
function inn_j(oj) result(jj)
  integer,intent(in) :: oj
  integer :: jj
  jj=oj+nj1-1
end function inn_j
function inn_k(ok) result(kk)
  integer,intent(in) :: ok
  integer :: kk
  kk=ok+nk1-1
end function inn_k

function out_i(ii) result(oi)
  integer,intent(in) :: ii
  integer :: oi
  oi=ii-ni1+1
end function out_i
function out_j(jj) result(oj)
  integer,intent(in) :: jj
  integer :: oj
  oj=jj-nj1+1
end function out_j
function out_k(kk) result(ok)
  integer,intent(in) :: kk
  integer :: ok
  ok=kk-nk1+1
end function out_k
function loct_i(n) result(i)
  integer,intent(in) :: n
  integer :: i
  i=n+nx1-1
end function loct_i
function loct_j(n) result(j)
  integer,intent(in) :: n
  integer :: j
  j=n+ny1-1
end function loct_j
function loct_k(n) result(k)
  integer,intent(in) :: n
  integer :: k
  k=n+nz1-1
end function loct_k

end module para_mod
