!******************************************************************************!
!*  This program locates the seismic source                                   *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

#include "mod_macdrp.h"
!-----------------------------------------------------------------------------
program seis3d_source
!-----------------------------------------------------------------------------

use constants_mod
use string_mod
use math_mod
use nfseis_mod
use para_mod
use grid_mod
use media_mod
use src_mod
use mpi_mod
implicit none

real(SP),allocatable ::        &
    fault_M0(:,:),             &
    dip(:),strike(:),rake(:,:)
real(SP) :: src_hyper_height
logical,allocatable ::         &
    fault_flag(:)
character (len=SEIS_STRLEN) :: filenm
integer i,j,k,n,npt,gi,gj,gk,n_i,n_j,n_k,si,sj,sk,n1,n2,n3
logical iflag
integer :: nloc(1)
real,allocatable :: gx(:),gy(:),gz(:)

!----------------------------------------------------------------------!

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)
call media_fnm_init(fnm_conf)
call src_fnm_init(fnm_conf)

! load grid
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

! source
call read_src_para(fnm_src_conf)

! find src index
do n=1,npt_fault
   nloc=minloc(abs(fault_loct(1,n)-gx))
   fault_indx(1,n)=loct_i(nloc(1))
   nloc=minloc(abs(fault_loct(2,n)-gy))
   fault_indx(2,n)=loct_j(nloc(1))
   nloc=minloc(abs(fault_loct(3,n)-gz))
   fault_indx(3,n)=loct_k(nloc(1))
   gi=fault_indx(1,n);gj=fault_indx(2,n)
   if (fault_loct(3,n)>src_hyper_height) then
      fault_indx(3,n)=n3-LenFD
   end if
   gk=fault_indx(3,n)

   if (     gi<ni1 .or. gi>n1-LenFD  &
       .or. gj<nj1 .or. gj>n2-LenFD  &
       .or. gk<nk1 .or. gk>n3-LenFD ) then
      print *, 'source out the model'
      print *, 'source point is ',fault_loct(:,n)
      print *, 'x scale of model:',gx(ni1),gx(n1-LenFD)
      print *, 'y scale of model:',gy(nj1),gy(n2-LenFD)
      print *, 'z scale of model:',gz(nk1),gz(n3-LenFD)
      stop 1
   end if
   fault_flag(n)=.false.
end do

! check
iflag=.false.
do n=1,npt_fault
   if ( fault_flag(n)) then
      iflag=.true.
      print *, 'n=',n,'loct=',fault_loct(:,n),'indx=',fault_indx(:,n)
   end if
end do
if (iflag) then
   print *, 'there are some source points left'
   stop 1
end if

! distribute to thread
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   filenm=get_fnm_srcnode(n_i,n_j,n_k)
   call nfseis_srcnode_create(filenm,ntwin,2*LenFD+1)
   call nfseis_srcnode_attput(filenm,                 &
        flag_stf_type, flag_mech_type,stf_t0, stf_alpha)
   call nfseis_varput(filenm,'stf_A',stf_A,           &
        (/1/),(/ntwin/),(/1/) )
   call nfseis_varput(filenm,'twin',twin,             &
        (/1,1/),(/2,ntwin/),(/1,1/) )

   n1=0
   do n=1,npt_fault
      gi=fault_indx(1,n);gj=fault_indx(2,n);gk=fault_indx(3,n)
   if (      gi>=ngx1 .and. gi<=ngx2                          &
       .and. gj>=ngy1 .and. gj<=ngy2                          &
       .and. gk>=ngz1 .and. gk<=ngz2 ) then
      n1=n1+1
      call nfseis_varput(filenm,'MomTxx',MomTxx(:,n),         &
           (/1,n1/),(/ntwin,1/),(/1,1/) )
      call nfseis_varput(filenm,'MomTyy',MomTyy(:,n),         &
           (/1,n1/),(/ntwin,1/),(/1,1/) )
      call nfseis_varput(filenm,'MomTzz',MomTzz(:,n),         &
           (/1,n1/),(/ntwin,1/),(/1,1/) )
      call nfseis_varput(filenm,'MomTxy',MomTxy(:,n),         &
           (/1,n1/),(/ntwin,1/),(/1,1/) )
      call nfseis_varput(filenm,'MomTxz',MomTxz(:,n),         &
           (/1,n1/),(/ntwin,1/),(/1,1/) )
      call nfseis_varput(filenm,'MomTyz',MomTyz(:,n),         &
           (/1,n1/),(/ntwin,1/),(/1,1/) )
      call nfseis_varput(filenm,'RuptT',RuptT(n),             &
           (/n1/),(/1/),(/1/) )
      call nfseis_varput(filenm,'fault_indx',                 &
           (/ out_i(swmpi_locli(gi,n_i)),                     &
              out_j(swmpi_loclj(gj,n_j)),                     &
              out_k(swmpi_loclk(gk,n_k)) /) ,                 &
           (/1,n1/),(/SEIS_GEO,1/),(/1,1/) )
      call nfseis_varput(filenm,'fault_loct',fault_loct(:,n), &
           (/1,n1/),(/SEIS_GEO,1/),(/1,1/) )
   end if
   end do
end do
end do
end do

call src_destroy
call dealloc_local

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

!*************************************************************************
!*                 PART-I  src alloc and dealloc                         *
!*************************************************************************

subroutine alloc_local
  allocate(fault_flag(npt_fault)); fault_flag=.true.
  allocate(fault_M0(ntwin,npt_fault)); fault_M0=0.0_SP
  allocate(dip(npt_fault)); dip=0.0_SP
  allocate(strike(npt_fault)); strike=0.0_SP
  allocate(rake(ntwin,npt_fault)); rake=0.0_SP
end subroutine alloc_local

subroutine dealloc_local
  if (allocated(fault_loct)) deallocate(fault_loct)
  if (allocated(fault_flag)) deallocate(fault_flag)
  if (allocated(fault_M0)) deallocate(fault_M0)
  if (allocated(dip)) deallocate(dip)
  if (allocated(strike)) deallocate(strike)
  if (allocated(rake)) deallocate(rake)
end subroutine dealloc_local

!*************************************************************************
!*                    PART-II  src init and distrib                      *
!*************************************************************************

subroutine read_src_para(fnm_conf)
character (len=*) :: fnm_conf
character (len=SEIS_STRLEN) :: srcMech,str,stf_type
real(SP) Mxx,Myy,Mzz,Mxy,Mxz,Myz,Mtm
real(SP) strike_pi,dip_pi,rake_pi,d2m
integer fid,i,n,m
fid=1001
open(fid,file=trim(fnm_conf),status="old")

call string_conf(fid,1,'distance2meter',2,d2m)

call string_conf(fid,1,"ntwin",2,ntwin)
call string_conf(fid,1,"npt_fault",2,npt_fault)
call string_conf(fid,1,"src_mech",2,srcMech)

call alloc_local
call src_alloc(npt_fault,ntwin)

do n=1,ntwin
   call string_conf(fid,1,"twin",2*n,twin(1,n))
   call string_conf(fid,1,"twin",2*n+1,twin(2,n))
end do

call string_conf(fid,1,'src_hyper_height',2,src_hyper_height)
src_hyper_height=src_hyper_height*d2m
call string_conf(fid,1,"<Epicenter",2,str)
do n=1,npt_fault
   read(fid,*) fault_loct(:,n), fault_M0(:,n), RuptT(n), &
               strike(n),dip(n),rake(:,n)
   fault_loct(1:2,n)=fault_loct(1:2,n)*PI/180.0_SP
   fault_loct(3,n)=fault_loct(3,n)*d2m
end do

if (trim(srcMech)=="double") then
   !fault_M0=fault_M0*1.0e-6  ! N m = kg m/s^2 m -> kg km/s^2 km
   flag_mech_type=SIG_MECH_DOUBLE
elseif (trim(srcMech)=="explosive") then
   !fault_M0=fault_M0*1.0e-6  ! N m = kg m/s^2 m -> kg km/s^2 km
   flag_mech_type=SIG_MECH_EXPLOSIVE
elseif (trim(srcMech)=="force") then
   !fault_M0=fault_M0*1.0e-3  ! N = kg m/s^2 -> kg km/s^2
   flag_mech_type=SIG_MECH_FORCE
end if

do n=1,npt_fault
   dip_pi=dip(n)/180.0_SP*PI
   strike_pi=strike(n)/180.0_SP*PI
do m=1,ntwin
   if (flag_mech_type==SIG_MECH_DOUBLE) then
      rake_pi=rake(m,n)/180.0_SP*PI
      Mxx=-(sin(dip_pi)*cos(rake_pi)*sin(2.0*strike_pi)   &
               +sin(2.0*dip_pi)*sin(rake_pi)*sin(strike_pi)**2)
      Myy=sin(dip_pi)*cos(rake_pi)*sin(2.0*strike_pi)     &
               -sin(2.0*dip_pi)*sin(rake_pi)*cos(strike_pi)**2
      !Mzz=sin(2.0*dip_pi)*sin(rake_pi)
      Mzz=-(Mxx+Myy)
      Mxy=sin(dip_pi)*cos(rake_pi)*cos(2.0*strike_pi)     &
               +0.5*sin(2.0*dip_pi)*sin(rake_pi)*sin(2.0*strike_pi)
      Mxz=-(cos(dip_pi)*cos(rake_pi)*cos(strike_pi)       &
               +cos(2.0*dip_pi)*sin(rake_pi)*sin(strike_pi))
      Myz=-(cos(dip_pi)*cos(rake_pi)*sin(strike_pi)       &
               -cos(2.0*dip_pi)*sin(rake_pi)*cos(strike_pi))
      !Mxz=-Mxz;Mxy=-Mxy !for upward positive z axis
      Mtm=Mxx; Mxx=Myy; Myy=Mtm
      Mtm=Mxz; Mxz=-Myz; Myz=-Mtm
   elseif (flag_mech_type==SIG_MECH_EXPLOSIVE) then
      Mxx=1.0_SP/sqrt(3.0_SP); Myy=1.0_SP/sqrt(3.0_SP); Mzz=1.0_SP/sqrt(3.0_SP)
      Mxy=0.0_SP; Mxz=0.0_SP; Myz=0.0_SP
   elseif (flag_mech_type==SIG_MECH_FORCE) then
      Mxx=strike(n)
      Myy=dip(n)
      Mzz=rake(m,n)
      Mxy=0.0_SP; Mxz=0.0_SP; Myz=0.0_SP
   end if

   MomTxx(m,n)=Mxx*fault_M0(m,n)
   MomTyy(m,n)=Myy*fault_M0(m,n)
   MomTzz(m,n)=Mzz*fault_M0(m,n)
   MomTxy(m,n)=Mxy*fault_M0(m,n)
   MomTxz(m,n)=Mxz*fault_M0(m,n)
   MomTyz(m,n)=Myz*fault_M0(m,n)
end do ! ntwin
end do ! npt_fault

call string_conf(fid,1,"stf",2,stf_type)
select case (trim(stf_type))
case ('bell_int')
     stf_t0=twin(2,1)-twin(1,1)
     stf_alpha=0.0_SP
     flag_stf_type=SIG_SVF_BELL
case ('gauss_int')
     call string_conf(fid,1,"stf_para",2,stf_alpha)
     call string_conf(fid,1,"stf_para",3,stf_t0)
     flag_stf_type=SIG_SVF_GAUSS
case ('triangle_int')
     stf_t0=twin(2,1)-twin(1,1)
     stf_alpha=0.0_SP
     flag_stf_type=SIG_SVF_TRIANGLE
case ('ricker')
     call string_conf(fid,1,"stf_para",2,stf_alpha)
     call string_conf(fid,1,"stf_para",3,stf_t0)
     flag_stf_type=SIG_STF_RICKER
case ('B(shift)')
     call string_conf(fid,1,"stf_t0",2,stf_t0)
     flag_stf_type=SIG_STF_BSHIFT
case ('gauss')
     call string_conf(fid,1,"stf_para",2,stf_alpha)
     call string_conf(fid,1,"stf_para",3,stf_t0)
     flag_stf_type=SIG_STF_GAUSS
case ('bell')
     stf_t0=twin(2,1)-twin(1,1)
     stf_alpha=0.0_SP
     flag_stf_type=SIG_STF_BELL
case default
     print *, "Have you told me how to generate the STF of ", trim(stf_type)
     stop 1
end select

do i=1,ntwin
   stf_A(i)=1.0
   !stf_A(i)=cal_stf_norm(twin(2,i)-twin(1,i),stept)
end do

close(fid)
end subroutine read_src_para

end program seis3d_source
