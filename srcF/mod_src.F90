!******************************************************************************!
!*  This module is used for incorporating seismic source                      *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

#include "mod_macdrp.h"
!-----------------------------------------------------------------------------
module src_mod
!-----------------------------------------------------------------------------

use constants_mod
use string_mod
use math_mod
use para_mod
use mpi_mod
use grid_mod
use media_mod
use nfseis_mod
implicit none
private

public ::          &
  src_fnm_init,    &
  src_alloc,       &
  src_destroy,     &
  src_import,      &
  get_fnm_srcnode, &
  stf_t0,stf_alpha,&
  src_stress,      &
  src_force
public ::        &
  fun_bell,      &
  fun_ricker,    &
  fun_gauss

DEFFDWET
HOCWETL
HOCWETR
!----------------------------------------------------------------------!

integer,parameter,public ::  &
     SIG_SVF_BELL     =100,  &
     SIG_SVF_TRIANGLE =110,  &
     SIG_SVF_GAUSS    =120,  &
     SIG_STF_RICKER   =200,  &
     SIG_STF_BSHIFT   =210,  &
     SIG_STF_GAUSS    =220,  &
     SIG_STF_BELL     =230
integer,parameter,public ::  &
     SIG_MECH_FORCE     =10, &
     SIG_MECH_DOUBLE    =20, &
     SIG_MECH_EXPLOSIVE =30

integer,public ::                             &
     npt_fault,                               &
     ntwin
integer,allocatable,public ::                 &
    fault_indx(:,:)
real(SP),allocatable,public ::                &
    fault_loct(:,:),                          &
    RuptT(:),                                 &
    twin(:,:),                                &
    stf_A(:)
real(SP),dimension(:,:),allocatable,public :: &
     MomTxx,MomTyy,MomTzz,                    &
     MomTxy,MomTxz,MomTyz
character (len=SEIS_STRLEN),public ::         &
     fnm_src_conf,                            &
     pnm_src
integer,public ::                             &
     flag_stf_type,                           &
     flag_mech_type
real(SP) :: stf_t0,stf_alpha
real(SP) :: stf_disp
real(SP) :: gs_s0,gs_t0
#ifdef VERBOSE
integer fid_out
#endif

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!
subroutine src_fnm_init(filenm)
character (len=*) :: filenm
integer fid
fid=1001
open(fid,file=trim(filenm),status="old")
  call string_conf(fid,1,'fnm_src_conf',2,fnm_src_conf)
  call string_conf(fid,1,'pnm_src',2,pnm_src)
close(fid)
end subroutine src_fnm_init

!*************************************************************************
!*                 PART-I  src alloc and dealloc                         *
!*************************************************************************
subroutine src_alloc(npt_fault,ntwin)
  integer npt_fault,ntwin
  allocate(fault_indx(SEIS_GEO,npt_fault)); fault_indx=0
  allocate(fault_loct(SEIS_GEO,npt_fault)); fault_loct=0.0
  allocate(RuptT(npt_fault)); RuptT=0.0
  allocate(twin(2,ntwin))
  allocate(stf_A(ntwin)); stf_A=0.0
  allocate(MomTxx(ntwin,npt_fault)); MomTxx=0.0
  allocate(MomTyy(ntwin,npt_fault)); MomTyy=0.0
  allocate(MomTzz(ntwin,npt_fault)); MomTzz=0.0
  allocate(MomTxy(ntwin,npt_fault)); MomTxy=0.0
  allocate(MomTxz(ntwin,npt_fault)); MomTxz=0.0
  allocate(MomTyz(ntwin,npt_fault)); MomTyz=0.0
end subroutine src_alloc

subroutine src_destroy
  if (allocated(MomTxx)) deallocate(MomTxx)
  if (allocated(MomTyy)) deallocate(MomTyy)
  if (allocated(MomTzz)) deallocate(MomTzz)
  if (allocated(MomTxy)) deallocate(MomTxy)
  if (allocated(MomTxz)) deallocate(MomTxz)
  if (allocated(MomTyz)) deallocate(MomTyz)
#ifdef VERBOSE
  close(fid_out)
#endif
end subroutine src_destroy

!*************************************************************************
!*                    PART-II source io                                  *
!*************************************************************************

subroutine src_import
  integer n,i,j,k
  character (len=SEIS_STRLEN) :: filenm
  filenm=get_fnm_srcnode(thisid(1),thisid(2),thisid(3))
  call nfseis_diminfo(filenm,'npt',npt_fault)
  call nfseis_diminfo(filenm,'ntwin',ntwin)
  if (npt_fault<=0) return

  call src_alloc(npt_fault,ntwin)
  ! Mij
  call nfseis_varget(filenm,'MomTxx',MomTxx, &
       (/1,1/),(/ntwin,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'MomTyy',MomTyy, &
       (/1,1/),(/ntwin,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'MomTzz',MomTzz, &
       (/1,1/),(/ntwin,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'MomTxy',MomTxy, &
       (/1,1/),(/ntwin,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'MomTxz',MomTxz, &
       (/1,1/),(/ntwin,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'MomTyz',MomTyz, &
       (/1,1/),(/ntwin,npt_fault/),(/1,1/))
  !Other
  call nfseis_varget(filenm,'fault_indx',fault_indx, &
       (/1,1/),(/SEIS_GEO,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'fault_loct',fault_loct, &
       (/1,1/),(/SEIS_GEO,npt_fault/),(/1,1/))
  call nfseis_varget(filenm,'RuptT',RuptT, &
       (/1/),(/npt_fault/),(/1/))
  call nfseis_varget(filenm,'twin',twin,   &
       (/1,1/),(/2,ntwin/),(/1,1/))
  call nfseis_varget(filenm,'stf_A',stf_A, &
       (/1/),(/ntwin/),(/1/))
  !att
  call nfseis_attget(filenm,'flag_stf_type',flag_stf_type)
  call nfseis_attget(filenm,'flag_mech_type',flag_mech_type)
  call nfseis_attget(filenm,'stf_t0',stf_t0)
  call nfseis_attget(filenm,'stf_alpha',stf_alpha)
  !correct
  do n=1,npt_fault
     i=fault_indx(1,n); j=fault_indx(2,n); k=fault_indx(3,n)
     fault_indx(1,n)=inn_i(i)
     fault_indx(2,n)=inn_j(j)
     fault_indx(3,n)=inn_k(k)
  end do
  ! for gauss stf type
  gs_s0=0.0; gs_t0=0.0
#ifdef VERBOSE
  fid_out=9010
  open(fid_out,                                                                      &
       file='log_stf_'//trim(set_mpi_subfix(thisid(1),thisid(2),thisid(3)))//'.dat', &
       status='unknown')
  write(fid_out,*) 'stf=',flag_stf_type
  write(fid_out,*) 'mech=',flag_mech_type
  write(fid_out,*) 't0 and alpha=', stf_t0, stf_alpha
  write(fid_out,*) 'ntime  time  stf'
#endif
end subroutine src_import

!*************************************************************************
!*                    PART-II source time function                       *
!*************************************************************************
!--- source signal generator --
!    partial derivative with respect to time
function fun_gauss(t,a,t0) result(f)
  real(SP) :: t,t0,a,f
  if (t<=0.0 .or. t>=2.0*t0) then
     f=0.0
  else
     f=exp(-(t-t0)**2/(a**2))/(sqrt(PI)*a)
  end if
end function fun_gauss

function fun_ricker(t,fc,t0) result (v)
    real(SP) :: t,t0,fc,v
    real(SP) :: u,f0
    if (t<=0.0) then
       v=0.0; return
    end if
    f0=sqrt(PI)/2.0
    u=(t-t0)*2.0*PI*fc
    v=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
    !v=u*(1.5-u**2/4.0)*exp(-u**2.0/4.0)*f0*PI*fc
end function fun_ricker
function fun_ricker_deriv(t,fc,t0) result (v)
    real(SP) :: t,t0,fc,v
    real(SP) :: u,f0
    if (t<=0.0) then
       v=0.0; return
    end if
    f0=sqrt(PI)/2.0
    u=(t-t0)*2.0*PI*fc
    !f=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
    v=u*(1.5-u**2/4.0)*exp(-u**2.0/4.0)*f0*PI*fc
end function fun_ricker_deriv
function fun_bell(t,riset) result(v)
  real(SP) :: t,riset,v
  if (t>0.0 .and. t<riset) then
     v=(1.0 - cos(2*PI*t/riset))/riset
  else
     v=0.0
  end if
end function fun_bell
function fun_bell_deriv(t,riset) result(v)
  real t,riset,v
  if (t>0.0 .and. t<riset) then
     v=2.0*PI*sin(2*PI*t/riset)/riset
  else
     v=0.0
  end if
end function fun_bell_deriv
function fun_bell_int(t,riset) result(v)
  real(SP) :: t,riset,v
  if (t<=0.0) then
     v=0.0
  elseif (t<riset) then
     v=t/riset - sin(2*PI*t/riset)/(2.0*pi)
  else
     v=1.0
  end if
end function fun_bell_int
function fun_triangle(t,riset) result(v)
  real(SP) :: t,t0,riset,v
  t0=riset/2.0;
  if (t>riset) then
     v=0.0
  elseif (t>t0) then
     v=2.0/t0-t/(t0**2)
  elseif (t>0.0) then
     v=t/(t0**2)
  else
     v=0.0
  end if
end function fun_triangle
function fun_triangle_int(t,riset) result(v)
  real(SP) :: t,t0,riset,v
  t0=riset/2.0;
  if (t>riset) then
     v=1.0
  elseif (t>t0) then
     v=-0.5*t**2/(t0**2)+2.0*t/t0-1.0
  elseif (t>0.0) then
     v=0.5*t**2/(t0**2)
  else
     v=0.0
  end if
end function fun_triangle_int
function fun_bshift(t,riset,t0) result(v)
    real(SP) :: t,t0,riset,bshift,v
    bshift=riset/2.0;
    if (t>0.0 .and. t<riset) then
       v=0.5/t0/(cosh((t-bshift)/t0)**2.0)
    else
       v=0.0
    end if
end function fun_bshift


function cal_svf_rate(t,riset,a0) result(rate)
real(SP) :: t,riset,a0,rate
select case(flag_stf_type)
case (SIG_SVF_BELL)
   rate= fun_bell(t,riset)/a0 
case (SIG_STF_BELL)
   rate= fun_bell_deriv(t,riset)/a0 
case (SIG_SVF_TRIANGLE)
   rate= fun_triangle(t,riset)/a0 
case (SIG_SVF_GAUSS)
   rate=fun_gauss(t,stf_alpha,stf_t0)/a0
case (SIG_STF_BSHIFT)
   rate= fun_bshift(t,riset,stf_t0)/a0 
case (SIG_STF_RICKER)
   rate=fun_ricker_deriv(t,stf_alpha,stf_t0)/a0
case default
   print *, "Have you told me how to generate the SVF of ", flag_stf_type
   stop 1
end select
if (abs(rate)<SEIS_ZERO) rate=0.0
end function cal_svf_rate
function cal_stf_disp(t,riset,a0) result(rate)
integer,parameter :: gs_nt=10
real(SP) :: t,riset,a0,rate
real(SP) :: dt,s1,s2
integer n
select case(flag_stf_type)
case (SIG_SVF_BELL)
   rate= fun_bell_int(t,riset)/a0 
case (SIG_STF_BELL)
   rate= fun_bell(t,riset)/a0 
case (SIG_STF_RICKER)
   rate=fun_ricker(t,stf_alpha,stf_t0)/a0
case (SIG_STF_GAUSS)
   rate=fun_gauss(t,stf_alpha,stf_t0)/a0
case (SIG_SVF_GAUSS)
   if (abs(t-gs_t0)<SEIS_ZERO) then
      rate=gs_s0
   else
      dt=(t-gs_t0)/gs_nt
      s1=fun_gauss(gs_t0,stf_alpha,stf_t0)
      s2=fun_gauss(t,stf_alpha,stf_t0)
      gs_s0=gs_s0+0.5*dt*(s1+s2)
      do n=1,gs_nt-1
         s1=fun_gauss(gs_t0+n*dt,stf_alpha,stf_t0)/a0
         gs_s0=gs_s0+dt*s1
      end do
      rate=gs_s0
      gs_t0=t
   end if
case default
   print *, "Have you told me how to generate the STF of ", flag_stf_type
   stop 1
end select
if (abs(rate)<SEIS_ZERO) rate=0.0
end function cal_stf_disp

subroutine cal_norm_delt(delt,x0,y0,z0,x,y,z,rx0,ry0,rz0)
  real,dimension(-LenFD:LenFD,-LenFD:LenFD,-LenFD:LenFD) :: delt
  real,dimension(-LenFD:LenFD) :: x,y,z
  real x0,y0,z0,rx0,ry0,rz0
  real D1,D2,D3
  integer i,j,k
  delt=0.0
  do k=-LenFD,LenFD
  do j=-LenFD,LenFD
  do i=-LenFD,LenFD
     D1=fun_guass(x(i)-x0,rx0)
     D2=fun_guass(y(j)-y0,ry0)
     D3=fun_guass(z(k)-z0,rz0)
     delt(i,j,k)=D1*D2*D3
  end do
  end do
  end do
  delt=delt/sum(delt)
end subroutine cal_norm_delt
function fun_guass(r,r0) result(f)
  real r,r0,f
  f=exp(-r**2/(r0**2))
end function fun_guass

!*************************************************************************
!*                    PART-I  src coupled                                *
!*************************************************************************
#ifdef SrcSmooth
subroutine src_stress(Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,tinc,stept)
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer ntime
real(SP) :: tinc,stept
integer i,j,k,n,m,si,sj,sk
real(SP) :: t,rate,Mxx,Mxy,Mxz,Myy,Myz,Mzz,Vr
real(SP),dimension(-LenFD:LenFD,-LenFD:LenFD,-LenFD:LenFD) :: normd
real x0,y0,z0,V

if (      flag_mech_type/=SIG_MECH_DOUBLE  &
    .and. flag_mech_type/=SIG_MECH_EXPLOSIVE ) return

t=(ntime+tinc)*stept

do n=1,npt_fault
   x0=fault_loct(1,n); y0=fault_loct(2,n); z0=fault_loct(3,n)
   i=fault_indx(1,n); j=fault_indx(2,n); k=fault_indx(3,n)
call cal_norm_delt(normd, x0,y0,z0,                          &
     x(i-LenFD:i+LenFD),y(j-LenFD:j+LenFD),z(k-LenFD:k+LenFD), &
     (x(i+1)-x(i-1))/2.0,(y(j+1)-y(j-1))/2.0,(z(k+1)-z(k-1))/2.0 )
if (freenode .and. k+LenFD>nk2)  then
   normd=normd/sum(normd(:,:,-LenFD:nk2-k))
end if

do m=1,ntwin
#ifdef SrcTensorHook
   rate=cal_svf_rate(t-RuptT(n)-twin(1,m), &
                  (twin(2,m)-twin(1,m)), stf_A(m) )
#endif
#ifdef SrcTensorMomentum
   rate=cal_stf_disp(t-RuptT(n)-twin(1,m), &
                  (twin(2,m)-twin(1,m)), stf_A(m) )
#endif
#ifdef VERBOSE
if (n==1) then
   write(fid_out,'(i5,2es12.5)') ntime, t, rate
end if
#endif
   Mxx=MomTxx(m,n); Myy=MomTyy(m,n); Mzz=MomTzz(m,n)
   Mxy=MomTxy(m,n); Mxz=MomTxz(m,n); Myz=MomTyz(m,n)
do sk=-LenFD,LenFD
do sj=-LenFD,LenFD
do si=-LenFD,LenFD
   i=fault_indx(1,n)+si; j=fault_indx(2,n)+sj; k=fault_indx(3,n)+sk
   if (      i<=ni2 .and. i>=ni1  &
       .and. j<=nj2 .and. j>=nj1  &
       .and. k<=nk2 .and. k>=nk1 ) then
       V=stepx*stepy*stepz*z(k)**2*xsin(i)
       Vr=rate/V
       Txx(i,j,k)=Txx(i,j,k)-Mxx*normd(si,sj,sk)*Vr
       Tyy(i,j,k)=Tyy(i,j,k)-Myy*normd(si,sj,sk)*Vr
       Tzz(i,j,k)=Tzz(i,j,k)-Mzz*normd(si,sj,sk)*Vr
       Txy(i,j,k)=Txy(i,j,k)-Mxy*normd(si,sj,sk)*Vr
       Txz(i,j,k)=Txz(i,j,k)-Mxz*normd(si,sj,sk)*Vr
       Tyz(i,j,k)=Tyz(i,j,k)-Myz*normd(si,sj,sk)*Vr
   end if
end do
end do
end do  ! loop_of_sindx

end do ! ntwin
end do ! loop_of_npt_fault
end subroutine src_stress
#else
subroutine src_stress(Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,tinc,stept)
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer ntime
real(SP) :: tinc,stept
integer i,j,k,n,m,si,sj,sk
real(SP) :: t,rate,Mxx,Mxy,Mxz,Myy,Myz,Mzz,Vr
real x0,y0,z0,V

if (      flag_mech_type/=SIG_MECH_DOUBLE  &
    .and. flag_mech_type/=SIG_MECH_EXPLOSIVE ) return

t=(ntime+tinc)*stept

do n=1,npt_fault
   x0=fault_loct(1,n); y0=fault_loct(2,n); z0=fault_loct(3,n)
   i=fault_indx(1,n); j=fault_indx(2,n); k=fault_indx(3,n)

do m=1,ntwin
#ifdef SrcTensorHook
   rate=cal_svf_rate(t-RuptT(n)-twin(1,m), &
                  (twin(2,m)-twin(1,m)), stf_A(m) )
#endif
#ifdef SrcTensorMomentum
   rate=cal_stf_disp(t-RuptT(n)-twin(1,m), &
                  (twin(2,m)-twin(1,m)), stf_A(m) )
#endif
#ifdef VERBOSE
if (n==1) then
   write(fid_out,'(i5,2es12.5)') ntime, t, rate
end if
#endif
   Mxx=MomTxx(m,n); Myy=MomTyy(m,n); Mzz=MomTzz(m,n)
   Mxy=MomTxy(m,n); Mxz=MomTxz(m,n); Myz=MomTyz(m,n)
   if (      i<=ni2 .and. i>=ni1  &
       .and. j<=nj2 .and. j>=nj1  &
       .and. k<=nk2 .and. k>=nk1 ) then
       V=stepx*stepy*stepz*z(k)**2*xsin(i)
       Vr=rate/V
       Txx(i,j,k)=Txx(i,j,k)-Mxx*Vr
       Tyy(i,j,k)=Tyy(i,j,k)-Myy*Vr
       Tzz(i,j,k)=Tzz(i,j,k)-Mzz*Vr
       Txy(i,j,k)=Txy(i,j,k)-Mxy*Vr
       Txz(i,j,k)=Txz(i,j,k)-Mxz*Vr
       Tyz(i,j,k)=Tyz(i,j,k)-Myz*Vr
   end if
end do ! ntwin
end do ! loop_of_npt_fault
end subroutine src_stress
#endif

subroutine src_force(Vx,Vy,Vz,ntime,tinc,stept)
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Vx,Vy,Vz
integer ntime
real(SP) :: tinc,stept

real(SP),dimension(-2*LenFD:2*LenFD,-2*LenFD:2*LenFD,-2*LenFD:2*LenFD) :: &
     normd
integer :: n,m,i,j,k,si,sj,sk,mi,mj,mk
real(SP) :: Mxx,Myy,Mzz
real(SP) :: t,rate,fx,fy,fz
real(SP) :: rrho
real(SP) x0,y0,z0,V

if (flag_mech_type/=SIG_MECH_FORCE) return

t=(ntime+tinc)*stept

do n=1,npt_fault
   x0=fault_loct(1,n); y0=fault_loct(2,n); z0=fault_loct(3,n)
   i=fault_indx(1,n); j=fault_indx(2,n); k=fault_indx(3,n)
#ifdef SrcSmooth
call cal_norm_delt(normd, x0,y0,z0,                          &
     x(i-LenFD:i+LenFD),y(j-LenFD:j+LenFD),z(k-LenFD:k+LenFD), &
     (x(i+1)-x(i-1))/2.0,(y(j+1)-y(j-1))/2.0,(z(k+1)-z(k-1))/2.0 )
if (freenode .and. k+LenFD>nk2)  then
   normd=normd/sum(normd(:,:,-LenFD:nk2-k))
end if
#else
   normd=0.0_SP; normd(0,0,0)=1.0_SP
#endif
   si=fault_indx(1,n); sj=fault_indx(2,n); sk=fault_indx(3,n)
do m=1,ntwin
   rate=cal_stf_disp(t-RuptT(n)-twin(1,m), &
                  (twin(2,m)-twin(1,m)), stf_A(m) )
#ifdef VERBOSE
if (n==1) then
   write(fid_out,'(i5,2es12.5)') ntime, t, rate
end if
#endif
   Mxx=MomTxx(m,n); Myy=MomTyy(m,n); Mzz=MomTzz(m,n)
do mk=-LenFD,LenFD
do mj=-LenFD,LenFD
do mi=-LenFD,LenFD
   i=mi+si; j=mj+sj; k=mk+sk
if (      i<=ni2 .and. i>=ni1  &
    .and. j<=nj2 .and. j>=nj1  &
    .and. k<=nk2 .and. k>=nk1 ) then
   !if (freenode .and. k==nk2) cycle
   V=stepx*stepy*stepz*z(k)**2*xsin(i)
   rrho=rate/(V*rho(i,j,k))
   fx=Mxx*rrho*normd(mi,mj,mk)
   fy=Myy*rrho*normd(mi,mj,mk)
   fz=Mzz*rrho*normd(mi,mj,mk)
   if (freenode .and. k==nk2) then
      fx=2.0_SP*fx;fy=2.0_SP*fy;fz=2.0_SP*fz
   end if
   Vx(i,j,k)=Vx(i,j,k)+fx
   Vy(i,j,k)=Vy(i,j,k)+fy
   Vz(i,j,k)=Vz(i,j,k)+fz
end if !i[jk]
end do !mi
end do !mj
end do !mk

end do ! ntwin
end do ! loop_of_npt_fault
!call src_charac(t)
end subroutine src_force

subroutine src_charac(t)
use macdrp_mod, only : &
    TxSrc,TySrc,TzSrc
real(SP) :: t

real(SP),dimension(SEIS_GEO) :: e3
real(SP) :: disp,Mxx,Myy,Mzz,area,fct
integer :: n,m,si,sj,sk,nf

if (.not. freenode) return
if (flag_mech_type/=SIG_MECH_FORCE) return

m=1 ! for ntwin

do nf=1,npt_fault
   si=fault_indx(1,nf); sj=fault_indx(2,nf); sk=fault_indx(3,nf)
   if (sk/=nk2) cycle
   disp=cal_stf_disp(t-RuptT(nf)-twin(1,m), &
                  (twin(2,m)-twin(1,m)), stf_A(m) )
   area=z(sk)**2*xsin(si)*stepx*stepy*stepz
   Mxx=MomTxx(m,nf); Myy=MomTyy(m,nf); Mzz=MomTzz(m,nf)
   TxSrc(si,sj)=disp*Mxx/area
   TySrc(si,sj)=disp*Myy/area
   TzSrc(si,sj)=disp*Mzz/area
end do ! loop_of_npt_fault
end subroutine src_charac

function get_fnm_srcnode(n_i,n_j,n_k) result(filenm)
  use io_mod, only : io_out_pattern
  integer n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_src)                          &
       //'/'//'src'                             &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'.nc'
end function get_fnm_srcnode

end module src_mod

