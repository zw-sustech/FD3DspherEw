!******************************************************************************!
!*  This program inits the media of the 3D geology model                      *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
program seis3d_media
!-----------------------------------------------------------------------------

#include "mod_macdrp.h"

use constants_mod
use string_mod, only : string_conf
use math_mod
use para_mod
use mpi_mod
use nfseis_mod
use grid_mod
use media_mod
use custom_mod
use io_mod

implicit none

integer :: NGRIDX,NGRIDY,NGRIDZ

type STRUCT_MEDIALAYERED
     logical :: yes
     integer :: ni,nj,nk
     real(SP),dimension(:),pointer :: x,y
     real(SP),dimension(:,:,:),pointer :: z
     real(SP),dimension(:,:,:),pointer :: h
     real(SP),dimension(:,:),pointer :: Vp,Vs,rho
     real(SP),dimension(:,:),pointer :: Qs
     real(SP) :: Qf0
end type STRUCT_MEDIALAYERED

type STRUCT_MEDIAPOLY
     logical :: yes
     character (len=SEIS_STRLEN) :: fnm
     integer :: imax,jmax,kmax,nod
     integer :: ni,nj,nd
     real(SP),dimension(:),pointer :: gx,gy
     real(SP),dimension(:),pointer :: x,y
     real(SP),dimension(:,:,:),pointer :: d
     real(SP),dimension(:,:,:,:),pointer :: aVp,aVs,aRho
     real(SP),dimension(:,:,:,:),pointer :: aQs
     real(SP) :: Qf0
end type STRUCT_MEDIAPOLY

type (STRUCT_MEDIALAYERED) :: L
type (STRUCT_MEDIAPOLY) :: W,P

integer,dimension(SEIS_GEO) ::         &
     subs,subc,subt,                   &
     subs_x1,subc_x1, subs_x2,subc_x2, &
     subs_y1,subc_y1, subs_y2,subc_y2, &
     subs_z1,subc_z1, subs_z2,subc_z2
integer,dimension(LenFD) :: indx_x1,indx_x2,indx_y1,indx_y2,indx_z1,indx_z2
integer :: i,j,k,n_i,n_j,n_k
character (len=SEIS_STRLEN) :: filenm

real(DP),dimension(:,:,:),allocatable :: gx,gy,gz
real(SP) :: ztopo
real(SP) :: Vp,dtmax
real(DP) :: dtlocal
integer,dimension(SEIS_GEO) :: dtnode,dtindx
real(SP) :: dtLe,dtmaxVp,dtmaxL

!-----------------------------------------------------------------------------

!fnm_media_conf='SeisMedia.conf'

call get_conf_name(fnm_conf)
call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)
call media_fnm_init(fnm_conf)
call media_alloc
call grid_alloc

!-----------------------------------------------------------------------------

print *, 'init media ...'
call init_media(fnm_media_conf)

!-----------------------------------------------------------------------------
print *, 'calculate effective media ...'
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=dims(3)-1,0,-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call grid_import
   if (n_k==dims(3)-1) then
      ztopo=z(nk2)
   end if

  if (P%yes) call mediapoly_read(P)

  if (W%yes) call mediapoly_read(W)

  call effmedia_eval

  call media_extend
  filenm=swmpi_rename_fnm(pnm_media,fnm_media)
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  call media_skel(filenm)
  call nfseis_varput(filenm,'rho',rho,subs,subc,subt)
  call nfseis_varput(filenm,'mu',mu,subs,subc,subt)
  call nfseis_varput(filenm,'lambda',lambda,subs,subc,subt)
#ifdef WITHQS
  call nfseis_varput(filenm,'Qs',Qs,subs,subc,subt)
#endif

end do
end do
end do

!-----------------------------------------------------------------------------
! exchange
subs_x1=(/ ni2+1,ny1  ,nz1  /); subc_x1=(/ LenFD , ny    , nz    /)
subs_x2=(/ nx1  ,ny1  ,nz1  /); subc_x2=(/ LenFD , ny    , nz    /)
subs_y1=(/ nx1  ,nj2+1,nz1  /); subc_y1=(/ nx    , LenFD , nz    /)
subs_y2=(/ nx1  ,ny1  ,nz1  /); subc_y2=(/ nx    , LenFD , nz    /)
subs_z1=(/ nx1  ,ny1  ,nk2+1/); subc_z1=(/ nx    , ny    , LenFD /)
subs_z2=(/ nx1  ,ny1  ,nz1  /); subc_z2=(/ nx    , ny    , LenFD /)
indx_x1=(/ (i,i=ni1,ni1+LenFD-1) /)
indx_x2=(/ (i,i=ni2-LenFD+1,ni2) /)
indx_y1=(/ (j,j=nj1,nj1+LenFD-1) /)
indx_y2=(/ (j,j=nj2-LenFD+1,nj2) /)
indx_z1=(/ (k,k=nk1,nk1+LenFD-1) /)
indx_z2=(/ (k,k=nk2-LenFD+1,nk2) /)
subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)

dtmax=1.0e10
allocate(gx(nx1:nx2,ny1:ny2,nz1:nz2)); gx=0.0_SP
allocate(gy(nx1:nx2,ny1:ny2,nz1:nz2)); gy=0.0_SP
allocate(gz(nx1:nx2,ny1:ny2,nz1:nz2)); gz=0.0_SP

write(*,*) "exchange media on boundary stencil ..."

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call media_import
! check
   call grid_import

   do k=nz1,nz2
   do j=ny1,ny2
   do i=nx1,nx2
      !gx(i,j,k)=z(k)*cos(real(y(j),DP))*cos(real(x(i),DP))
      !gy(i,j,k)=z(k)*cos(real(y(j),DP))*sin(real(x(i),DP))
      !gz(i,j,k)=z(k)*sin(real(y(j),DP))
      gx(i,j,k)=z(k)*cos(real(x(i)-PI/2,DP))*cos(real(y(j),DP))
      gy(i,j,k)=z(k)*cos(real(x(i)-PI/2,DP))*sin(real(y(j),DP))
      gz(i,j,k)=z(k)*sin(real(x(i)-PI/2,DP))
   end do
   end do
   end do
   
   do k=nk1,nk2
   do j=nj1,nj2
   do i=ni1,ni2
      Vp=sqrt( (lambda(i,j,k)+2.0*mu(i,j,k))/rho(i,j,k) )
      dtLe=min( dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k-1),gy(i  ,j  ,k-1),gz(i  ,j  ,k-1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j-1,k  ),gy(i  ,j-1,k  ),gz(i  ,j-1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i-1,j  ,k  ),gy(i-1,j  ,k  ),gz(i-1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) , &
                dist_point2plane(gx(i,j,k),gy(i,j,k),gz(i,j,k),             &
                   (/ gx(i+1,j  ,k  ),gy(i+1,j  ,k  ),gz(i+1,j  ,k  ) /),   &
                   (/ gx(i  ,j+1,k  ),gy(i  ,j+1,k  ),gz(i  ,j+1,k  ) /),   &
                   (/ gx(i  ,j  ,k+1),gy(i  ,j  ,k+1),gz(i  ,j  ,k+1) /)) )
      dtlocal=1.3/Vp * dtLe 
      if (dtlocal<dtmax) then
         dtmax=dtlocal
         dtnode=(/ n_i,n_j,n_k /)
         dtindx=(/ i,j,k /)
         dtmaxVp=Vp
         dtmaxL=dtLe
      end if
   end do
   end do
   end do
   ! to x1
if (n_i>0) then
   call swmpi_change_fnm(n_i-1,n_j,n_k)
   filenm=swmpi_rename_fnm(pnm_media,fnm_media)
   call nfseis_varput(filenm,'rho',rho(indx_x1,:,:),               &
        subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'mu'      ,mu(indx_x1,:,:)      ,     &
        subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'lambda'  ,lambda(indx_x1,:,:)  ,     &
        subs_x1,subc_x1,subt)
end if
! to x2
if (n_i<dims(1)-1) then
   call swmpi_change_fnm(n_i+1,n_j,n_k)
   filenm=swmpi_rename_fnm(pnm_media,fnm_media)
   call nfseis_varput(filenm,'rho', rho(indx_x2,:,:),              &
        subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'mu', mu(indx_x2,:,:),                &
        subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'lambda', lambda(indx_x2,:,:),        &
        subs_x2,subc_x2,subt)
end if
! to y1
if (n_j>0) then
   call swmpi_change_fnm(n_i,n_j-1,n_k)
   filenm=swmpi_rename_fnm(pnm_media,fnm_media)
   call nfseis_varput(filenm,'rho', rho(:,indx_y1,:),              &
        subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'mu', mu(:,indx_y1,:),                &
        subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'lambda', lambda(:,indx_y1,:),        &
        subs_y1,subc_y1,subt)
end if
! to y2
if (n_j<dims(2)-1) then
   call swmpi_change_fnm(n_i,n_j+1,n_k)
   filenm=swmpi_rename_fnm(pnm_media,fnm_media)
   call nfseis_varput(filenm,'rho', rho(:,indx_y2,:),              &
        subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'mu', mu(:,indx_y2,:),                &
        subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'lambda', lambda(:,indx_y2,:),        &
        subs_y2,subc_y2,subt)
end if
! to k1
if (n_k>0) then
   call swmpi_change_fnm(n_i,n_j,n_k-1)
   filenm=swmpi_rename_fnm(pnm_media,fnm_media)
   call nfseis_varput(filenm,'rho', rho(:,:,indx_z1),              &
        subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'mu', mu(:,:,indx_z1),                &
        subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'lambda', lambda(:,:,indx_z1),        &
        subs_z1,subc_z1,subt)
end if
! to k2
if (n_k<dims(3)-1) then
   call swmpi_change_fnm(n_i,n_j,n_k+1)
   filenm=swmpi_rename_fnm(pnm_media,fnm_media)
   call nfseis_varput(filenm,'rho', rho(:,:,indx_z2),              &
        subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'mu', mu(:,:,indx_z2),                &
        subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'lambda', lambda(:,:,indx_z2),        &
        subs_z2,subc_z2,subt)
end if
end do
end do
end do

print *, "Maximum allowed time step is", dtmax
write(*,"(a,3i5,a,3i5)") "located on", dtindx,' in thread', dtnode
print *, " Vp and dL are:", dtmaxVp,dtmaxL
if (dtmax<stept) then
   print *, "Serious Error: stept>dtmax", stept,dtmax
   stop 1
end if

call media_destroy
call grid_dealloc
deallocate(gx)
deallocate(gy)
deallocate(gz)

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

subroutine init_media(fnm_conf)
character (len=*),intent(in) :: fnm_conf
integer :: fid,imax,jmax,kmax,i,j,k
real(DP) :: d2m,v2ms,d2kgm3
character (len=SEIS_STRLEN) :: str,layer_type

fid=1001
open(fid,file=trim(fnm_conf),status="old")

call string_conf(fid,1,'half_sample_point',2,NGRIDX)
call string_conf(fid,1,'half_sample_point',3,NGRIDY)
call string_conf(fid,1,'half_sample_point',4,NGRIDZ)

! perturbed
P%yes=.false.
call string_conf(fid,1,'perturbation_nc',2,P%fnm)
if (trim(P%fnm)/='none') then
   P%yes=.true.
   call nfseis_diminfo(P%fnm,'theta',imax)
   call nfseis_diminfo(P%fnm,'phi',jmax)
   call nfseis_diminfo(P%fnm,'break',kmax)
   call nfseis_diminfo(P%fnm,'order',P%nod)
   allocate(P%gx(imax)); P%gx=0.0
   allocate(P%gy(jmax)); P%gy=0.0
   call nfseis_varget(P%fnm,'theta',P%gx,(/1/),(/imax/),(/1/))
   call nfseis_varget(P%fnm,'phi',P%gy,(/1/),(/jmax/),(/1/))
   P%gx=P%gx*PI/180.0_DP; P%gy=P%gy*PI/180.0_DP
   P%imax=imax; P%jmax=jmax; P%kmax=kmax
   imax=min(nx*max(2*NGRIDX,1),imax)
   jmax=min(ny*max(2*NGRIDY,1),jmax)
   kmax=min(nz*max(2*NGRIDZ,1),kmax)
   call mediapoly_alloc(P,imax,jmax,kmax)
end if

! ncinterp
W%yes=.false.
call string_conf(fid,1,'media_interp_nc',2,W%fnm)
if (trim(W%fnm)/='none') then
   W%yes=.true.
   call nfseis_diminfo(W%fnm,'theta',imax)
   call nfseis_diminfo(W%fnm,'phi',jmax)
   call nfseis_diminfo(W%fnm,'break',kmax)
   call nfseis_diminfo(W%fnm,'order',W%nod)
   allocate(W%gx(imax)); W%gx=0.0
   allocate(W%gy(jmax)); W%gy=0.0
   call nfseis_varget(W%fnm,'theta',W%gx,(/1/),(/imax/),(/1/))
   call nfseis_varget(W%fnm,'phi',W%gy,(/1/),(/jmax/),(/1/))
   W%gx=W%gx*PI/180.0_DP; W%gy=W%gy*PI/180.0_DP
   W%imax=imax; W%jmax=jmax; W%kmax=kmax
   imax=min(nx*max(2*NGRIDX,1),imax)
   jmax=min(ny*max(2*NGRIDY,1),jmax)
   kmax=min(nz*max(2*NGRIDZ,1),kmax)
   call mediapoly_alloc(W,imax,jmax,kmax)
   return
end if

! layered
L%yes=.false.
call string_conf(fid,1,'nlayer',2,kmax)
if (kmax<1) return

L%yes=.true.
call string_conf(fid,1,'distance2meter',2,d2m)
call string_conf(fid,1,'velocity2m/s',2,v2ms)
call string_conf(fid,1,'density2kg/m^3',2,d2kgm3)
!call string_conf(fid,1,'layer_sealevel',2,layer_sealevel)
call string_conf(fid,1,'imax',2,imax)
call string_conf(fid,1,'jmax',2,jmax)
L%ni=imax; L%nj=jmax; L%nk=kmax
allocate(L%x(imax))
allocate(L%y(jmax))
allocate(L%z(imax,jmax,kmax))
allocate(L%h(imax,jmax,kmax-1))
allocate(L%Vp (2,kmax))
allocate(L%Vs (2,kmax))
allocate(L%rho(2,kmax))
allocate(L%Qs (2,kmax))

do k=1,kmax
   call string_conf(fid,1,'vp',2*k   ,L%Vp(1,k))
   call string_conf(fid,1,'vp',2*k+1 ,L%Vp(2,k))
   call string_conf(fid,1,'vs',2*k   ,L%Vs(1,k))
   call string_conf(fid,1,'vs',2*k+1 ,L%Vs(2,k))
   call string_conf(fid,1,'rho',2*k  ,L%rho(1,k))
   call string_conf(fid,1,'rho',2*k+1,L%rho(2,k))
#ifdef WITHQS
   call string_conf(fid,1,'Q',2*k  ,  L%Qs(1,k))
   call string_conf(fid,1,'Q',2*k+1,  L%Qs(2,k))
#endif
end do
#ifdef WITHQS
call string_conf(fid,1,'Qf0',2, L%Qf0)
#endif

L%Vp(1,1)=L%Vp(2,1); L%Vs(1,1)=L%Vs(2,1)
L%rho(1,1)=L%rho(2,1); L%Qs(1,1)=L%Qs(2,1)
L%Vp=L%Vp*v2ms; L%Vs=L%Vs*v2ms; L%rho=L%rho*d2kgm3

! check
if ( any(L%rho<=0) .or. any(L%Vp<=0) .or. any(L%Vs<0) ) then
   print *, 'media parameter is negtive'
   print *, 'Vp=',L%Vp
   print *, 'Vs=',L%Vs
   print *, 'rho=',L%rho
   stop 1
end if
if ( any(L%Vp**2<=2.0*L%Vs**2) ) then
   print *, 'Vp**2 < 2Vs**2'
   print *, 'Vp=',L%Vp
   print *, 'Vs=',L%Vs
   stop 1
end if

call string_conf(fid,1,'layer_meaning',2,layer_type)
call string_conf(fid,1,'<Layer',2,str)
do j=1,jmax
do i=1,imax     
    read(fid,*) L%x(i),L%y(j),( L%z(i,j,k), k=1,kmax )
end do 
end do

select case (trim(layer_type))
case ('axis')
     do k=1,kmax-1
        L%h(:,:,k)=L%z(:,:,k)-L%z(:,:,k+1)
     end do
case ('depth')
     do k=1,kmax-1
        L%h(:,:,k)=L%z(:,:,k+1)-L%z(:,:,k)
     end do
case ('thickness')
     do k=1,kmax-1
        L%h(:,:,k)=L%z(:,:,k)
     end do
case default
     print *, "layer_meaning can only take axis, depth or thickness"
     stop 1
end select

do k=1,kmax-1
   if (any(L%h(:,:,k)<0.0)) then
      print *, 'thickness of layer',k,' should be not less than 0'
      print *, minloc(L%h(:,:,k))
      stop 1
   end if
end do
L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP
!L%z=L%z+layer_sealevel;
L%h=L%h*d2m

close(fid)

if (W%yes .and. L%yes) then
   print *, "cann't deal with poly and layered media at same time currently"
   stop 1
end if
end subroutine init_media

subroutine mediapoly_alloc(P,imax,jmax,kmax)
type(STRUCT_MEDIAPOLY) :: P
integer,intent(in) :: imax,jmax,kmax

if (associated(P%d)) then
if (size(P%d,1)<imax .or. size(P%d,2)<jmax .or. size(P%d,3)>=kmax) then
   deallocate(P%d); deallocate(P%aVp); deallocate(P%aVs); deallocate(P%aRho)
end if
end if

P%ni=imax;P%nj=jmax;P%nd=kmax
if (.not. associated(P%d)) then
   allocate(P%x(imax)); P%x=0.0
   allocate(P%y(jmax)); P%y=0.0
   allocate(P%d(kmax,imax,jmax)); P%d=0.0
   allocate(P%aVp(P%nod,kmax,imax,jmax)); P%aVp=0.0
   allocate(P%aVs(P%nod,kmax,imax,jmax)); P%aVs=0.0
   allocate(P%aRho(P%nod,kmax,imax,jmax)); P%aRho=0.0
end if
end subroutine mediapoly_alloc

subroutine mediapoly_read(P)
type(STRUCT_MEDIAPOLY) :: P
integer :: i1,i2,j1,j2,k1,k2,imax,jmax,kmax
integer :: indx(1)
indx=maxloc(P%gx,P%gx<=max(x(nx1),P%gx(1)));i1=max(indx(1),1)
indx=minloc(P%gx,P%gx>=min(x(nx2),P%gx(P%imax)));i2=indx(1)
indx=maxloc(P%gy,P%gy<=max(y(ny1),P%gy(1)));j1=max(indx(1),1)
indx=minloc(P%gy,P%gy>=min(y(ny2),P%gy(P%jmax)));j2=indx(1)
k1=1; k2=P%kmax
imax=i2-i1+1; jmax=j2-j1+1; kmax=k2-k1+1
call mediapoly_alloc(P,imax,jmax,kmax)

P%x(1:imax)=P%gx(i1:i2); P%y(1:jmax)=P%gy(j1:j2)
call nfseis_varget( P%fnm,'depth',P%d(1:kmax,1:imax,1:jmax), &
     (/k1,i1,j1/),(/kmax,imax,jmax/),(/1,1,1/) )
call nfseis_varget( P%fnm,'coefs_of_Vp',P%aVp(:,1:kmax,1:imax,1:jmax), &
     (/1,k1,i1,j1/),(/P%nod,kmax,imax,jmax/),(/1,1,1,1/) )
call nfseis_varget( P%fnm,'coefs_of_Vs',P%aVs(:,1:kmax,1:imax,1:jmax), &
     (/1,k1,i1,j1/),(/P%nod,kmax,imax,jmax/),(/1,1,1,1/) )
call nfseis_varget( P%fnm,'coefs_of_rho',P%aRho(:,1:kmax,1:imax,1:jmax), &
     (/1,k1,i1,j1/),(/P%nod,kmax,imax,jmax/),(/1,1,1,1/) )
end subroutine mediapoly_read

subroutine dealloc_media
if (L%yes) then
   deallocate(L%x)
   deallocate(L%y)
   deallocate(L%h)
   deallocate(L%Vp)
   deallocate(L%Vs)
   deallocate(L%rho)
end if
end subroutine dealloc_media

subroutine effmedia_eval

real(SP),dimension(-NGRIDX:NGRIDX) :: xvec
real(SP),dimension(-NGRIDY:NGRIDY) :: yvec
real(SP),dimension(-NGRIDZ:NGRIDZ) :: zvec

real(SP),dimension(2) :: lam, miu, Qatt, Vp,Vs,dens
real(SP) :: rho0,mu0,lam0
real(SP) :: x0,y0,z0,d0
real(SP) :: dx1,dx2,dy1,dy2,dz1,dz2
integer :: nsamp
logical :: flag_water
integer :: i,j,k,mi,mj,mk,n,m
integer :: indx(1)

!layer
integer :: i1,i2,j1,j2,n1,n2
real(SP) :: h,d
real(SP) :: L1,L2
!poly
integer :: wi,wj
!perturbed
real(SP) :: dVp,dVs,dRho

!imax=L%ni; jmax=L%nj; kmax=L%nk
nsamp=max(2*NGRIDX,1)*max(2*NGRIDY,1)*max(2*NGRIDZ,1)

do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni2

   rho0=0.0;mu0=0.0;lam0=0.0
   flag_water=.false.

   xvec(0)=x(i); yvec(0)=y(j); zvec(0)=z(k)
   if (NGRIDX>0) then
      dx1=(x(i)-x(i-1))/2.0/NGRIDX; dx2=(x(i+1)-x(i))/2.0/NGRIDX
      xvec(-NGRIDX:-1)=(/-NGRIDX:-1/)*dx1+0.5*dx1+x(i);
      xvec(1:NGRIDX)=(/1:NGRIDX/)*dx2-0.5*dx2+x(i)
   end if
   if (NGRIDY>0) then
      dy1=(y(j)-y(j-1))/2.0/NGRIDY; dy2=(y(j+1)-y(j))/2.0/NGRIDY
      yvec(-NGRIDY:-1)=(/-NGRIDY:-1/)*dy1+0.5*dy1+y(j);
      yvec(1:NGRIDY)=(/1:NGRIDY/)*dy2-0.5*dy2+y(k)
   end if
   if (NGRIDZ>0) then
      dz1=(z(k)-z(k-1))/2.0/NGRIDZ; dz2=(z(k+1)-z(k))/2.0/NGRIDZ
      zvec(-NGRIDZ:-1)=(/-NGRIDZ:-1/)*dz1+0.5*dz1+z(k);
      zvec(1:NGRIDZ)=(/1:NGRIDZ/)*dz2-0.5*dz2+z(k)
   end if

do mk=-NGRIDZ,NGRIDZ
   if (NGRIDZ/=0 .and. mk==0) cycle
do mj=-NGRIDY,NGRIDY
   if (NGRIDY/=0 .and. mj==0) cycle
do mi=-NGRIDX,NGRIDX
   if (NGRIDX/=0 .and. mi==0) cycle
   x0=xvec(mi); y0=yvec(mj); z0=zvec(mk);
   d0=max(ztopo-z0,0.0);
   !x0=min(max(x0,L%x(1)),L%x(imax))
   !y0=min(max(y0,L%y(1)),L%y(jmax))
! background model
if (W%yes) then
   indx=minloc(abs(W%x(1:W%ni)-x0)); wi=indx(1)
   indx=minloc(abs(W%y(1:W%nj)-y0)); wj=indx(1)
   do n=1,W%nd
      if (abs(W%d(n,wi,wj)-d0)<=SEIS_ZERO) then
         if (n==1) then
            Vp(:)=ppval(W%aVp(:,n,wi,wj),d0-W%d(n,wi,wj))
            Vs(:)=ppval(W%aVs(:,n,wi,wj),d0-W%d(n,wi,wj))
            dens(:)=ppval(W%arho(:,n,wi,wj),d0-W%d(n,wi,wj))
         else
            Vp(1)  =ppval(W%aVp (:,n-1,wi,wj),d0-W%d(n-1,wi,wj));
            Vs(1)  =ppval(W%aVs (:,n-1,wi,wj),d0-W%d(n-1,wi,wj));
            dens(1)=ppval(W%arho(:,n-1,wi,wj),d0-W%d(n-1,wi,wj));
            Vp(2)  =ppval(W%aVp (:,n,wi,wj),d0-W%d(n,wi,wj));
            Vs(2)  =ppval(W%aVs (:,n,wi,wj),d0-W%d(n,wi,wj));
            dens(2)=ppval(W%arho(:,n,wi,wj),d0-W%d(n,wi,wj));
         end if
         exit
      elseif (W%d(n,wi,wj)>d0) then
         Vp(:)  =ppval(W%aVp (:,n-1,wi,wj),d0-W%d(n-1,wi,wj))
         Vs(:)  =ppval(W%aVs (:,n-1,wi,wj),d0-W%d(n-1,wi,wj))
         dens(:)=ppval(W%arho(:,n-1,wi,wj),d0-W%d(n-1,wi,wj))
         exit
      elseif (n==W%nd) then
         Vp(:)  =ppval(W%aVp (:,n,wi,wj),d0-W%d(n,wi,wj))
         Vs(:)  =ppval(W%aVs (:,n,wi,wj),d0-W%d(n,wi,wj))
         dens(:)=ppval(W%arho(:,n,wi,wj),d0-W%d(n,wi,wj))
      end if
   end do
end if
if (L%yes) then
   d=0.0
   call indx_locate_1d(x0,L%x,i1,i2)
   call indx_locate_1d(y0,L%y,j1,j2)
   do n=1,L%nk-1
      h=interp_2d(L%x(i1:i2),L%y(j1:j2),L%h(i1:i2,j1:j2,n),2,2,x0,y0)
      if (h<=SEIS_ZERO) then
         if (n==L%nk-1) then
            Vp(:)=L%Vp(2,L%nk)
            Vs(:)=L%Vs(2,L%nk)
            dens(:)=L%rho(2,L%nk)
         else
            cycle
         end if
      end if
      d=d+h
      if (abs(d-d0)<=SEIS_ZERO) then
         n1=n+1;
         do m=n+1,L%nk-1
            h=interp_2d(L%x(i1:i2),L%y(j1:j2),L%h(i1:i2,j1:j2,n),2,2,x0,y0)
            if (h>SEIS_ZERO) then
                n2=m
                exit
            end if
            if (m==L%nk-1) n2=m+1
         end do
         Vp(1)=L%Vp(1,n1); Vp(2)=L%Vp(2,n2)
         Vs(1)=L%Vs(1,n1); Vs(2)=L%Vs(2,n2)
         dens(1)=L%rho(1,n1); dens(2)=L%rho(2,n2)
         exit
      elseif (d>d0) then
         L1=(d-d0)/h; L2=1.0-L1
         Vp(:)=L%Vp(2,n)*L1+L%Vp(1,n+1)*L2
         Vs(:)=L%Vs(2,n)*L1+L%Vs(1,n+1)*L2
         dens(:)=L%rho(2,n)*L1+L%rho(1,n+1)*L2
         exit
      end if
   end do
end if

! perturbed
if (P%yes) then
   indx=minloc(abs(P%x(1:P%ni)-x0)); wi=indx(1)
   indx=minloc(abs(P%y(1:P%nj)-y0)); wj=indx(1)
   do n=2,P%nd
      if (P%d(n,wi,wj)>d0) then
         dVp =ppval(P%aVp (:,n-1,wi,wj),d0-P%d(n-1,wi,wj))
         dVs =ppval(P%aVs (:,n-1,wi,wj),d0-P%d(n-1,wi,wj))
         dRho=ppval(P%arho(:,n-1,wi,wj),d0-P%d(n-1,wi,wj))
         exit
      elseif (n==P%nd) then
         dVp =ppval(P%aVp (:,n,wi,wj),d0-P%d(n,wi,wj))
         dVs =ppval(P%aVs (:,n,wi,wj),d0-P%d(n,wi,wj))
         dRho=ppval(P%arho(:,n,wi,wj),d0-P%d(n,wi,wj))
      end if
   end do
   Vp=Vp*(1.0+dVp)
   Vs=Vs*(1.0+dVs)
   dens=dens*(1.0+dRho)
end if

   ! to lame parameters
   miu=dens*Vs*Vs
   lam=Vp*Vp*dens - 2.0*miu

   !accumulate
   rho0=rho0+0.5*(dens(1)+dens(2))
   if (miu(1)<=SEIS_ZERO .or. miu(2)<=SEIS_ZERO) then
      flag_water=.true.
   else
      mu0 = mu0+0.5*(1.0/miu(1)+1.0/miu(2))
   end if
   !mu0 = mu0+0.5*(1.0/miu(1)+1.0/miu(2))
   lam0=lam0+0.5*(1.0/lam(1)+1.0/lam(2))

end do !mi
end do !mj
end do !mk

    rho(i,j,k)=rho0/nsamp
    if (flag_water) then
       mu(i,j,k)=0.0
    else
       mu(i,j,k)=nsamp/mu0
    end if
    lambda(i,j,k)=nsamp/lam0

end do !i
end do !j
end do !k
end subroutine effmedia_eval

function ppval(pp,x0) result(f)
real(SP),dimension(:),intent(in) :: pp
real(SP),intent(in) :: x0
real(SP) :: f
f=pp(1)*x0**3+pp(2)*x0**2+pp(3)*x0+pp(4)
end function ppval

subroutine media_extend
  call extend_equal(rho)
  call extend_equal(mu)
  call extend_equal(lambda)
#ifdef WITHQS
  call extend_equal(Qs)
#endif
end subroutine media_extend
subroutine extend_equal(w)
   real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: w
   integer i,j,k,n
   ! x1, x2
   do k=nz1,nz2
   do j=ny1,ny2
      do n=1,LenFD
         w(ni1-n,j,k)=w(ni1,j,k)
         w(ni2+n,j,k)=w(ni2,j,k)
      end do
   end do
   end do
   ! y1, y2
   do k=nz1,nz2
   do i=nx1,nx2
      do n=1,LenFD
         w(i,nj1-n,k)=w(i,nj1,k)
         w(i,nj2+n,k)=w(i,nj2,k)
      end do
   end do
   end do
   ! z1, z2
   do j=ny1,ny2
   do i=nx1,nx2
      do n=1,LenFD
         w(i,j,nk1-n)=w(i,j,nk1)
         w(i,j,nk2+n)=w(i,j,nk2)
      end do
   end do
   end do
end subroutine extend_equal
function dist_point2plane(x0,y0,z0,A,B,C) result (L)
real(DP),intent(in) :: x0,y0,z0
real(DP),dimension(SEIS_GEO),intent(in) :: A,B,C
real(DP) L
real(DP),dimension(SEIS_GEO) :: AB,AC,p
real(DP) c1,c2,c3,d

AB=B-A; AC=C-A
call times_product(AB,AC,p)
c1=p(1);c2=p(2);c3=p(3);
d=-dot_product(p,A)
L=abs( (c1*x0+c2*y0+c3*z0+d)/sqrt(dot_product(p,p)) )
end function dist_point2plane

subroutine media_skel(filenm)
character (len=*),intent(in) :: filenm
integer,dimension(SEIS_GEO) ::subs,subc,subt
   call nfseis_data_create(filenm,nx,ny,nz, &
                       "media generated by seis3d_media" )
   subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
   call nfseis_data_attput(filenm,        &
        subs,subc,subt,                   &
        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &
        nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &
        ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &
        ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,    &
        (/ngi1,ngi2,ngj1,ngj2,ngk1,ngk2/) )
   call nfseis_data_addvar(filenm,'rho')
   call nfseis_data_addvar(filenm,'lambda')
   call nfseis_data_addvar(filenm,'mu')
#ifdef WITHQS
   call nfseis_data_addvar(filenm,'Qs')
#endif
end subroutine

subroutine indx_locate_1d(x0,x,i1,i2)
real(SP),intent(in) :: x0
real(SP),dimension(:),intent(in) :: x
integer,intent(out) :: i1,i2
!logical,intent(in) :: backward
integer :: indx(1)

!if (backward) then
!else
  indx=minloc(x,x>=x0)
  i2=max(indx(1),2)
  i1=i2-1
!end if
end subroutine indx_locate_1d

end program seis3d_media

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
