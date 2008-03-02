!******************************************************************************!
!*  This module is used for output data                                       *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
module io_mod
!-----------------------------------------------------------------------------

use constants_mod, only : SEIS_STRLEN,SEIS_GEO
use string_mod, only : string_conf
use para_mod
use nfseis_mod
use mpi_mod
implicit none

private
public :: io_init,io_pt_read,io_pt_import,io_snap_read,io_snap_locate
public :: io_rest_export,io_rest_import
public :: io_seismo_init,io_seismo_put,io_seismo_close
public :: io_wave_export,io_wave_close
public :: io_enum,io_out_pattern
public :: io_delete
public :: get_fnm_seismo,get_fnm_seismoinfo,       &
          get_fnm_snapnode,get_fnm_snapnode_n,     &
          get_fnm_snapcoord,get_fnm_snap,          &
          get_fnm_recv,get_fnm_line
public :: seismo_on_this,retrieve_recvline
public :: retrieve_snap_seis,retrieve_snap_time
public :: corr_subse,corr_indx
public :: read_nc_list

character (len=SEIS_STRLEN),public ::        &
    fnm_log

!---------------------------------------------
integer,public :: num_recv,num_line,num_snap,num_pt

integer,public :: pt_tinv
integer,allocatable,public :: &
     pt_indx(:,:),            &
     pt_id  (:,:),            &
     num_line_pt(:)
real(SP),allocatable,public ::    &
     pt_xyz(:,:)
real(SP),public :: topo_hyper_height

integer,allocatable,public :: &
     snap_id(:),              &
     snap_tinv(:),            &
     snap_tcnt(:),            &
     snap_gsubs(:,:),         &
     snap_gsubt(:,:),         &
     snap_gsube(:,:),         &
     snap_gsubc(:,:),         &
     snap_subs(:,:),          &
     snap_subt(:,:),          &
     snap_sube(:,:),          &
     snap_subc(:,:),          &
     snap_ncid(:),            &
     snap_vxid(:),            &
     snap_vyid(:),            &
     snap_vzid(:),            &
     snap_tid(:),             &
     snap_nid(:),             &
     norm_ncid(:),            &
     norm_vxid(:),            &
     norm_vyid(:),            &
     norm_vzid(:),            &
     norm_tid(:),             &
     norm_nid(:),             &
     Tsh_ncid(:),             &
     Tsh_vxid(:),             &
     Tsh_vyid(:),             &
     Tsh_vzid(:),             &
     Tsh_tid(:),              &
     Tsh_nid(:)
logical,allocatable,public :: &
     snap_oflag(:),           &
     stress_out(:)

character (len=SEIS_STRLEN),public ::           &
    pnm_out,pnm_seismo,pnm_seismoinfo,pnm_snap, &
    fnm_recv, prefix_line, prefix_snap,         &
    pnm_rest, fnm_rest, fnm_rest_point

!---------------------------------------------
integer pt_ncid,pt_tid,pt_vxid,pt_vyid,pt_vzid
integer rest_tinv,run_from_rest
integer :: nt_dyn_rest

!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
subroutine io_init(fnm_input)
character (len=*) :: fnm_input
integer fid

fid=1001
open(fid,file=trim(fnm_input),status="old")
!-- log file --
call string_conf(fid,1,'fnm_log',2,fnm_log)

!-- restart --
call string_conf(fid,1,'pnm_rest',2,pnm_rest)
call string_conf(fid,1,'rest_tinv',2,rest_tinv)
call string_conf(fid,1,'run_from_rest',2,run_from_rest)
call string_conf(fid,1,'fnm_rest_point',2,fnm_rest_point)
nt_dyn_rest=0

call string_conf(fid,1,'pnm_output',2,pnm_out)
call string_conf(fid,1,'pnm_seismo',2,pnm_seismo)
call string_conf(fid,1,'pnm_seismoinfo',2,pnm_seismoinfo)
call string_conf(fid,1,'pnm_snap',2,pnm_snap)
call string_conf(fid,1,'fnm_recv',2,fnm_recv)
call string_conf(fid,1,'prefix_line',2,prefix_line)
call string_conf(fid,1,'prefix_snap',2,prefix_snap)

!fnm_seismo='seis.nc'
!fnm_seismoinfo='seisinfo.nc'
!pref_snap='snap_'
!fnm_snapinfo='snapinfo.nc'
!fnm_recv='recv.nc'
!pref_line='line_'

close(fid)

!if (masternode) then
!   open(fid,file=trim(fnm_log),status="unknown")
!     write(fid,*) "# seis3d_wave run time log"
!   close(fid)
!end if
end subroutine io_init

!---------------------------------------------------------------------------

subroutine io_pt_read(fnm_input)
character (len=*) :: fnm_input
real(SP),dimension(SEIS_GEO) :: xyz0,dxyz
integer fid,n,m,npt

fid=1001
open(fid,file=trim(fnm_input),status="old")
call string_conf(fid,1,'number_of_recv',2,num_recv)
call string_conf(fid,1,'number_of_inline',2,num_line)
call string_conf(fid,1,'tinv_of_seismo',2,pt_tinv)
call string_conf(fid,1,'topo_hyper_height',2,topo_hyper_height)
num_pt=num_recv
do n=1,num_line
   call string_conf(fid,1,trim(io_enum('line_',n)),8,m)
   num_pt=num_pt+m
end do
call alloc_pt(num_pt,num_line)
!-- recv --
npt=0
do n=1,num_recv
   npt=npt+1
   do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('recv_',n)),m+1,pt_xyz(m,npt))
   end do
   pt_id(:,npt)=(/ n,0 /)
end do
!-- inline --
do n=1,num_line
   do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('line_',n)),m+1,xyz0(m))
      call string_conf(fid,1,trim(io_enum('line_',n)),m+1+SEIS_GEO,dxyz(m))
   end do
   call string_conf(fid,1,trim(io_enum('line_',n)),1+1+2*SEIS_GEO,num_line_pt(n))
   do m=1,num_line_pt(n)
      npt=npt+1
      pt_xyz(:,npt)=xyz0+dxyz*(m-1)
      pt_id(:,npt)=(/ m,n /)
   end do
end do
pt_xyz(1:2,:)=pt_xyz(1:2,:)*PI/180.0_SP
end subroutine io_pt_read
subroutine alloc_pt(npt,nline)
 integer npt,nline
 allocate(pt_xyz(SEIS_GEO,npt))
 allocate(pt_id (       2,npt))
 allocate(num_line_pt(nline))
end subroutine alloc_pt

subroutine io_pt_import
integer n
character (len=SEIS_STRLEN) :: filenm
filenm= get_fnm_seismoinfo(pnm_seismoinfo,thisid(1),thisid(2),thisid(3))
call nfseis_diminfo(trim(filenm),'num_pt',num_pt)
call nfseis_attget(trim(filenm),'tinv',pt_tinv)
if (num_pt>0) then
   allocate(pt_indx(SEIS_GEO,num_pt))
   call nfseis_varget(trim(filenm),'indx', pt_indx, &
        (/1,1/),(/SEIS_GEO,num_pt/),(/1,1/))
do n=1,num_pt
   pt_indx(1,n)=inn_i(pt_indx(1,n))
   pt_indx(2,n)=inn_j(pt_indx(2,n))
   pt_indx(3,n)=inn_k(pt_indx(3,n))
end do
end if
end subroutine io_pt_import

!---------------------------------------------------------------------------

subroutine io_snap_read(fnm_input)
character (len=*) :: fnm_input
integer fid,n,m

fid=1001
open(fid,file=trim(fnm_input),status="old")
!-- output --
call string_conf(fid,1,'number_of_snap',2,num_snap)
call alloc_snap(num_snap)
do n=1,num_snap
do m=1,SEIS_GEO
 call string_conf(fid,1,trim(io_enum('snap_',n)),m+1,snap_subs(m,n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),m+1+SEIS_GEO,snap_subc(m,n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),m+1+2*SEIS_GEO,snap_subt(m,n))
end do
 call string_conf(fid,1,trim(io_enum('snap_',n)),1+1+3*SEIS_GEO,snap_tinv(n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),2+1+3*SEIS_GEO,snap_tcnt(n))
 call string_conf(fid,1,trim(io_enum('snap_',n)),3+1+3*SEIS_GEO,stress_out(n))
 call corr_indx(snap_subs(:,n),snap_subc(:,n),snap_subt(:,n),snap_sube(:,n))
 snap_id(n)=n
end do
!snap_sube=snap_subs+snap_subt*(snap_subc-1)
close(fid)
end subroutine io_snap_read

subroutine alloc_snap(nsnap)
  integer nsnap
  allocate(snap_gsubs(SEIS_GEO,nsnap))
  allocate(snap_gsube(SEIS_GEO,nsnap))
  allocate(snap_gsubt(SEIS_GEO,nsnap))
  allocate(snap_gsubc(SEIS_GEO,nsnap))
  allocate(snap_subs(SEIS_GEO,nsnap))
  allocate(snap_sube(SEIS_GEO,nsnap))
  allocate(snap_subt(SEIS_GEO,nsnap))
  allocate(snap_subc(SEIS_GEO,nsnap))
  allocate(snap_tinv(nsnap))
  allocate(snap_tcnt(nsnap))
  allocate(snap_id(nsnap))
  allocate(snap_ncid(nsnap))
  allocate(snap_vxid(nsnap))
  allocate(snap_vyid(nsnap))
  allocate(snap_vzid(nsnap))
  allocate(snap_tid(nsnap))
  allocate(snap_nid(nsnap))
  allocate(snap_oflag(nsnap)); snap_oflag=.false.
  allocate(norm_ncid(nsnap))
  allocate(norm_vxid(nsnap))
  allocate(norm_vyid(nsnap))
  allocate(norm_vzid(nsnap))
  allocate(norm_tid(nsnap))
  allocate(norm_nid(nsnap))
  allocate(Tsh_ncid(nsnap))
  allocate(Tsh_vxid(nsnap))
  allocate(Tsh_vyid(nsnap))
  allocate(Tsh_vzid(nsnap))
  allocate(Tsh_tid(nsnap))
  allocate(Tsh_nid(nsnap))
  allocate(stress_out(nsnap)); stress_out=.false.
end subroutine alloc_snap

subroutine io_snap_locate
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
integer,dimension(SEIS_GEO) :: gsubs,gsube
integer n,nsnap,tinv,tcnt

nsnap=0
do n=1,num_snap
   subs=snap_subs(:,n);subc=snap_subc(:,n);subt=snap_subt(:,n)
   sube=snap_sube(:,n)
   tinv=snap_tinv(n)
   tinv=snap_tinv(n)
   tcnt=snap_tcnt(n)
   ! convert into this thread
if (      subs(1)<=ngi2 .and. sube(1)>=ngi1 &
    .and. subs(2)<=ngj2 .and. sube(2)>=ngj1 &
    .and. subs(3)<=ngk2 .and. sube(3)>=ngk1 ) then
   nsnap=nsnap+1
   call corr_subse(subs,subc,subt,sube)
   gsubs(1)=out_i(subs(1)); gsubs(2)=out_j(subs(2)); gsubs(3)=out_k(subs(3))
   gsube(1)=out_i(sube(1)); gsube(2)=out_j(sube(2)); gsube(3)=out_k(sube(3))
   snap_gsubs(:,nsnap)=gsubs
   snap_gsube(:,nsnap)=gsube
   snap_gsubc(:,nsnap)=subc
   snap_gsubt(:,nsnap)=subt
   snap_subs(1,nsnap)=swmpi_locli(subs(1),thisid(1))
   snap_subs(2,nsnap)=swmpi_loclj(subs(2),thisid(2))
   snap_subs(3,nsnap)=swmpi_loclk(subs(3),thisid(3))
   snap_sube(1,nsnap)=swmpi_locli(sube(1),thisid(1))
   snap_sube(2,nsnap)=swmpi_loclj(sube(2),thisid(2))
   snap_sube(3,nsnap)=swmpi_loclk(sube(3),thisid(3))
   snap_subc(:,nsnap)=subc
   snap_subt(:,nsnap)=subt
   snap_id(nsnap)=n
   snap_tinv(nsnap)=tinv
   snap_tcnt(nsnap)=tcnt
   stress_out(nsnap)=stress_out(n)
end if
end do
num_snap=nsnap
end subroutine io_snap_locate
!---------------------------------------------------------------------------

subroutine corr_indx(subs,subc,subt,sube)
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
   ! -1 global first index
   ! 0 index of source center
   ! -2 global index of free surface
   if (subs(1)==-1) then
       subs(1)=ni1
   else
       subs(1)=inn_i(subs(1))
   end if

   if (subs(2)==-1) then
       subs(2)=nj1
   else
       subs(2)=inn_j(subs(2))
   end if

   if (subs(3)==-2) then
      subs(3)=swmpi_globk(nk2,dims(3)-1)
   elseif (subs(3)==-1) then
      subs(3)=nk1
   else
      subs(3)=inn_k(subs(3))
   end if

   sube=subs+subt*(subc-1)
   if (subc(1)==-1) then
      sube(1)=swmpi_globi(ni2,dims(1)-1)
      subc(1)=(sube(1)-subs(1))/subt(1)+1
   end if
   if (subc(2)==-1) then
      sube(2)=swmpi_globj(nj2,dims(2)-1)
      subc(2)=(sube(2)-subs(2))/subt(2)+1
   end if
   if (subc(3)==-1) then
      sube(3)=swmpi_globk(nk2,dims(3)-1)
      subc(3)=(sube(3)-subs(3))/subt(3)+1
   end if
end subroutine corr_indx
subroutine corr_subse(subs,subc,subt,sube)
! in global mode
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
integer n
 if (ngi1>subs(1)) then
    n=ngi1-subs(1)
    if (mod(n,subt(1))==0) then
       subs(1)=ngi1
    else
       n=(n+subt(1)-1)/subt(1)
       subs(1)=n*subt(1)+subs(1)
    end if
 end if
 if (ngj1>subs(2)) then
    n=ngj1-subs(2)
    if (mod(n,subt(2))==0) then
       subs(2)=ngj1
    else
       n=(n+subt(2)-1)/subt(2)
       subs(2)=n*subt(2)+subs(2)
    end if
 end if
 if (ngk1>subs(3)) then
    n=ngk1-subs(3)
    if (mod(n,subt(3))==0) then
       subs(3)=ngk1
   else
       n=(n+subt(3)-1)/subt(3)
       subs(3)=n*subt(3)+subs(3)
   end if
 end if
 subc(1)=(ngi2-subs(1))/subt(1)
 subc(2)=(ngj2-subs(2))/subt(2)
 subc(3)=(ngk2-subs(3))/subt(3)
 sube(1)=min(sube(1),subc(1)*subt(1)+subs(1))
 sube(2)=min(sube(2),subc(2)*subt(2)+subs(2))
 sube(3)=min(sube(3),subc(3)*subt(3)+subs(3))
 subc=(sube-subs)/subt+1
end subroutine corr_subse
!---------------------------------------------------------------------------
subroutine io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
use mpi
integer ntime
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz
integer,dimension(SEIS_GEO) :: subs,subc,subt
character (len=SEIS_STRLEN) :: filenm
integer ierr,fid,n
  fid=5001
if (mod(ntime,10)==0) then
if (masternode) then
  open(fid,file=trim(fnm_rest_point),status='old')
  read(fid,*) nt_dyn_rest
  close(fid)
end if
  call MPI_BCAST(nt_dyn_rest,1,MPI_INTEGER,0,SWMPI_COMM,ierr)
end if

if (ntime/=nt_dyn_rest .and. mod(ntime,rest_tinv)/=0) return

  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  fnm_rest='rest_t'//trim(io_out_pattern(ntime,5))//'.nc'
  filenm=swmpi_rename_fnm(pnm_rest,fnm_rest)
  call nfseis_data_create( filenm,nx,ny,nz, &
       "Restart file generated by seis3d_wave" )
  call nfseis_data_addvar(filenm,'Vx')
  call nfseis_data_addvar(filenm,'Vy')
  call nfseis_data_addvar(filenm,'Vz')
  call nfseis_data_addvar(filenm,'Txx')
  call nfseis_data_addvar(filenm,'Tyy')
  call nfseis_data_addvar(filenm,'Tzz')
  call nfseis_data_addvar(filenm,'Txy')
  call nfseis_data_addvar(filenm,'Txz')
  call nfseis_data_addvar(filenm,'Tyz')
  call nfseis_varput( filenm,'Vx', Vx, subs,subc,subt)
  call nfseis_varput( filenm,'Vy', Vy, subs,subc,subt)
  call nfseis_varput( filenm,'Vz', Vz, subs,subc,subt)
  call nfseis_varput( filenm,'Txx',Txx,subs,subc,subt)
  call nfseis_varput( filenm,'Tyy',Tyy,subs,subc,subt)
  call nfseis_varput( filenm,'Tzz',Tzz,subs,subc,subt)
  call nfseis_varput( filenm,'Txy',Txy,subs,subc,subt)
  call nfseis_varput( filenm,'Txz',Txz,subs,subc,subt)
  call nfseis_varput( filenm,'Tyz',Tyz,subs,subc,subt)
  if (num_pt>0) ierr=nf90_sync(pt_ncid)
  do n=1,num_snap
     if (snap_oflag(n)) then
        ierr=nf90_sync(snap_ncid(n))
        if (stress_out(n)) then
           ierr=nf90_sync(norm_ncid(n))
           ierr=nf90_sync(Tsh_ncid(n))
        end if
     end if
  end do

  if (ntime-rest_tinv>0) then
     fnm_rest='rest_t'//trim(io_out_pattern(ntime-rest_tinv,5))//'.nc'
     filenm=swmpi_rename_fnm(pnm_rest,fnm_rest)
     call io_delete(filenm)
  end if
end subroutine io_rest_export
subroutine io_rest_import(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
integer ntime
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz
integer,dimension(SEIS_GEO) :: subs,subc,subt
character (len=SEIS_STRLEN) :: filenm
if (run_from_rest==0) then
  ntime=0
else
  ntime=run_from_rest
  subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  fnm_rest='rest_t'//trim(io_out_pattern(run_from_rest,5))//'.nc'
  filenm=swmpi_rename_fnm(pnm_rest,fnm_rest)
  call nfseis_varget( filenm,'Vx', Vx, subs,subc,subt)
  call nfseis_varget( filenm,'Vy', Vy, subs,subc,subt)
  call nfseis_varget( filenm,'Vz', Vz, subs,subc,subt)
  call nfseis_varget( filenm,'Txx',Txx,subs,subc,subt)
  call nfseis_varget( filenm,'Tyy',Tyy,subs,subc,subt)
  call nfseis_varget( filenm,'Tzz',Tzz,subs,subc,subt)
  call nfseis_varget( filenm,'Txy',Txy,subs,subc,subt)
  call nfseis_varget( filenm,'Txz',Txz,subs,subc,subt)
  call nfseis_varget( filenm,'Tyz',Tyz,subs,subc,subt)
end if
end subroutine io_rest_import
!---------------------------------------------------------------------------

subroutine io_seismo_init
if (num_pt<=0) return
if (run_from_rest==0) then
   call nfseis_seismo_init(                                             &
        get_fnm_seismo(pnm_out,thisid(1),thisid(2),thisid(3)),          &
        num_pt=num_pt,                                                  &
        ncid=pt_ncid,tid=pt_tid,vxid=pt_vxid,vyid=pt_vyid,vzid=pt_vzid)
else
   call nfseis_seismo_reinit(                                           &
        get_fnm_seismo(pnm_out,thisid(1),thisid(2),thisid(3)),          &
        ncid=pt_ncid,tid=pt_tid,vxid=pt_vxid,vyid=pt_vyid,vzid=pt_vzid)
end if
end subroutine io_seismo_init

subroutine io_seismo_put(Vx,Vy,Vz,ntime)
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Vx,Vy,Vz
integer ntime
real(SP) :: t
integer i,j,k,n,m,ierr

if (num_pt<=0) return

t=ntime*stept
if (mod(ntime,pt_tinv)==0) then
m=ntime/pt_tinv
do n=1,num_pt
   i=pt_indx(1,n);j=pt_indx(2,n);k=pt_indx(3,n)
   call nfseis_put(pt_ncid,pt_tid,t,(/m/),(/1/),(/1/))
   call nfseis_put(pt_ncid,pt_vxid,Vx(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_vyid,Vy(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
   call nfseis_put(pt_ncid,pt_vzid,Vz(i,j,k), &
        (/n,m/),(/1,1/),(/1,1/))
end do
end if

if (mod(ntime,pt_tinv*100)==0) ierr=nf90_sync(pt_ncid)
end subroutine io_seismo_put

subroutine io_seismo_close
if (num_pt<=0) return
call nfseis_close(pt_ncid)
end subroutine io_seismo_close

!---------------------------------------------------------------------------
subroutine io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Vx,Vy,Vz
real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
integer,dimension(SEIS_GEO) :: isubs,isube
!integer,dimension(SEIS_GEO) :: gsubs,gsubc,gsubt,gsube
integer ntime
real(SP) :: stept
character (len=SEIS_STRLEN) :: filenm
real(SP) :: t
integer n,n1
t=ntime*stept
do n=1,num_snap
if ( mod(ntime,snap_tinv(n))==0) then
   n1=mod(ntime/snap_tinv(n)-1,snap_tcnt(n))
   subs=snap_subs(:,n); subt=snap_subt(:,n); sube=snap_sube(:,n)
   subc=snap_subc(:,n)
   isubs(1) =out_i(subs(1));isubs(2)=out_j(subs(2));isubs(3)=out_k(subs(3))
   isube(1) =out_i(sube(1));isube(2)=out_j(sube(2));isube(3)=out_k(sube(3))

if ( n1==0 ) then
   filenm=get_fnm_snapnode(pnm_out,'snap_',n,ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_wave_init(filenm,'Vx','Vy','Vz',                &
          snap_ncid(n),snap_vxid(n),snap_vyid(n),snap_vzid(n), &
          snap_tid(n),snap_nid(n),stept*snap_tinv(n),          &
          snap_gsubs(:,n),                                     &
          snap_gsubc(:,n),                                     &
          snap_gsubt(:,n),                                     &
          snap_gsube(:,n),                                     &
                    isubs,                                     &
                     subc,                                     &
                     subt,                                     &
                    isube,                                     &
          "snap of wave feilds")
   snap_oflag(n)=.true.

   if (stress_out(n)) then
   filenm=get_fnm_snapnode(pnm_out,'norm_',n,                  &
          ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_wave_init(filenm,'Txx','Tyy','Tzz',             &
          norm_ncid(n),norm_vxid(n),norm_vyid(n),norm_vzid(n), &
          norm_tid(n),norm_nid(n),stept*snap_tinv(n),          &
          snap_gsubs(:,n),                                     &
          snap_gsubc(:,n),                                     &
          snap_gsubt(:,n),                                     &
          snap_gsube(:,n),                                     &
                    isubs,                                     &
                     subc,                                     &
                     subt,                                     &
                    isube,                                     &
          "snap of normal stress")
   filenm=get_fnm_snapnode(pnm_out,'shear_',n,             &
          ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_wave_init(filenm,'Txy','Txz','Tyz',         &
          Tsh_ncid(n),Tsh_vxid(n),Tsh_vyid(n),Tsh_vzid(n), &
          Tsh_tid(n),Tsh_nid(n),stept*snap_tinv(n),        &
          snap_gsubs(:,n),                                 &
          snap_gsubc(:,n),                                 &
          snap_gsubt(:,n),                                 &
          snap_gsube(:,n),                                 &
                    isubs,                                 &
                     subc,                                 &
                     subt,                                 &
                    isube,                                 &
          "snap of shear stress")
   end if
elseif ( run_from_rest>0 .and. ntime-run_from_rest<=snap_tinv(n) ) then
   filenm=get_fnm_snapnode(pnm_out,'snap_',n,                  &
          ntime,                                               &
          thisid(1),thisid(2),thisid(3))
   call nfseis_wave_reinit(filenm,'Vx','Vy','Vz',              &
          snap_ncid(n),snap_vxid(n),snap_vyid(n),snap_vzid(n), &
          snap_tid(n),snap_nid(n))
   snap_oflag(n)=.true.
   if (stress_out(n)) then
   filenm=get_fnm_snapnode(pnm_out,'norm_',n,                  &
          ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_wave_reinit(filenm,'Txx','Tyy','Tzz',           &
          norm_ncid(n),norm_vxid(n),norm_vyid(n),norm_vzid(n), &
          norm_tid(n),norm_nid(n))
   filenm=get_fnm_snapnode(pnm_out,'shear_',n,                 &
          ntime,thisid(1),thisid(2),thisid(3))
   call nfseis_wave_reinit(filenm,'Txy','Txz','Tyz',           &
          Tsh_ncid(n),Tsh_vxid(n),Tsh_vyid(n),Tsh_vzid(n),     &
          Tsh_tid(n),Tsh_nid(n))
   end if
end if

   call nfseis_put(snap_ncid(n),snap_tid(n),t,     &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(snap_ncid(n),snap_nid(n),ntime, &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(snap_ncid(n),snap_vxid(n),      &
        Vx(subs(1):sube(1):subt(1),                &
           subs(2):sube(2):subt(2),                &
           subs(3):sube(3):subt(3)),               &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(snap_ncid(n),snap_vyid(n),      &
        Vy(subs(1):sube(1):subt(1),                &
           subs(2):sube(2):subt(2),                &
           subs(3):sube(3):subt(3)),               &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(snap_ncid(n),snap_vzid(n),      &
        Vz(subs(1):sube(1):subt(1),                &
           subs(2):sube(2):subt(2),                &
           subs(3):sube(3):subt(3)),               &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   if (stress_out(n)) then
   call nfseis_put(norm_ncid(n),norm_tid(n),t,     &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(norm_ncid(n),norm_nid(n),ntime, &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(norm_ncid(n),norm_vxid(n),      &
        Txx(subs(1):sube(1):subt(1),               &
            subs(2):sube(2):subt(2),               &
            subs(3):sube(3):subt(3)),              &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(norm_ncid(n),norm_vyid(n),      &
        Tyy(subs(1):sube(1):subt(1),               &
            subs(2):sube(2):subt(2),               &
            subs(3):sube(3):subt(3)),              &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(norm_ncid(n),norm_vzid(n),      &
        Tzz(subs(1):sube(1):subt(1),               &
            subs(2):sube(2):subt(2),               &
            subs(3):sube(3):subt(3)),              &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(Tsh_ncid(n),Tsh_tid(n),t,     &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(Tsh_ncid(n),Tsh_nid(n),ntime, &
        (/n1+1/),(/1/),(/1/) )
   call nfseis_put(Tsh_ncid(n),Tsh_vxid(n),      &
        Txy(subs(1):sube(1):subt(1),             &
            subs(2):sube(2):subt(2),             &
            subs(3):sube(3):subt(3)),            &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(Tsh_ncid(n),Tsh_vyid(n),      &
        Txz(subs(1):sube(1):subt(1),             &
            subs(2):sube(2):subt(2),             &
            subs(3):sube(3):subt(3)),            &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   call nfseis_put(Tsh_ncid(n),Tsh_vzid(n),      &
        Tyz(subs(1):sube(1):subt(1),             &
            subs(2):sube(2):subt(2),             &
            subs(3):sube(3):subt(3)),            &
        (/1,1,1,n1+1/),(/subc,1/),(/1,1,1,1/) )
   end if

if ( n1==snap_tcnt(n)-1 ) then
   call nfseis_close(snap_ncid(n))
   if (stress_out(n)) then
       call nfseis_close(norm_ncid(n))
       call nfseis_close(Tsh_ncid(n))
   end if
   snap_oflag(n)=.false.
end if

end if ! output
end do
end subroutine io_wave_export

subroutine io_wave_close
integer n
do n=1,num_snap
if (snap_oflag(n)) &
   call nfseis_close(snap_ncid(n))
   if (stress_out(n)) then
       call nfseis_close(norm_ncid(n))
       call nfseis_close(Tsh_ncid(n))
   end if
end do
end subroutine io_wave_close

subroutine read_nc_list(pnm_list,fnm_list,varnm,var,subs1,subc1)
character (len=*) :: pnm_list,fnm_list,varnm
real(SP),dimension(:,:) :: var
integer,dimension(SEIS_GEO) :: subs1,subc1

character (len=SEIS_STRLEN) :: fnm_nc
integer,dimension(SEIS_GEO) :: subs,subc,subt
integer fid,ierr

integer i1,i2,j1,j2,k1,k2
integer gi1,gi2,gj1,gj2,gk1,gk2
integer x1,y1,z1

gi1=subs1(1); gi2=subs1(1)+subc1(1)-1
gj1=subs1(2); gj2=subs1(2)+subc1(2)-1
gk1=subs1(3); gk2=subs1(3)+subc1(3)-1

fid=4001
open(fid,file=trim(pnm_list)//'/'//trim(fnm_list),status='old')
do
   read(fid,'(a)',iostat=ierr) fnm_nc
   if (ierr<0) exit

   fnm_nc=trim(pnm_list)//'/'//trim(fnm_nc)
   call nfseis_attget(fnm_nc,'indx',subs)
   call nfseis_diminfo(fnm_nc,'I',i2)
   call nfseis_diminfo(fnm_nc,'J',j2)
   call nfseis_diminfo(fnm_nc,'K',k2)
   i1=subs(1);j1=subs(2);k1=subs(3)
   i2=i2+i1-1;j2=j2+j1-1;k2=k2+k1-1

if (      i1<=gi2 .and. i2>=gi1 &
    .and. j1<=gj2 .and. j2>=gj1 &
    .and. k1<=gk2 .and. k2>=gk1 ) then
   x1=max(gi1,i1); y1=max(gj1,j1); z1=max(gk1,k1)

   subs(1)=x1-i1+1;
   subs(2)=y1-j1+1;
   subs(3)=z1-k1+1;
   subc(1)=min(gi2,i2)-x1+1
   subc(2)=min(gj2,j2)-y1+1
   subc(3)=min(gk2,k2)-z1+1
   subt=1
   i1=x1-gi1+1; i2=i1+subc(1)-1
   j1=y1-gj1+1; j2=j1+subc(2)-1
   k1=z1-gk1+1; k2=k1+subc(3)-1
   call nfseis_varget(fnm_nc,trim(varnm),var(i1:i2,j1:j2), &
        subs,subc,subt)
end if

end do !read
close(fid)
end subroutine read_nc_list

function seismo_on_this(pnm_info,id,indx,n_i,n_j,n_k,npt) result(isIn)
character (len=*) :: pnm_info
integer id,indx,n_i,n_j,n_k,npt
logical isIn

integer n,num_pt
character (len=SEIS_STRLEN) :: filenm
integer,allocatable :: ids(:,:)

isIn=.false.
filenm= get_fnm_seismoinfo(pnm_info,n_i,n_j,n_k)
call nfseis_diminfo(trim(filenm),'num_pt',num_pt)
if (num_pt>0) then
   allocate(ids(2,num_pt))
   call nfseis_varget(trim(filenm),'id', ids, &
        (/1,1/),(/2,num_pt/),(/1,1/))
do n=1,num_pt
if (ids(2,n)==id .and. ids(1,n)==indx) then
   isIn=.true.
   npt=n
   exit
end if
end do
   deallocate(ids)
end if
end function seismo_on_this
subroutine retrieve_recvline(fnm_nc,id,Vx,Vy,Vz,nts,ntc,ntt,T)
character (len=*) :: fnm_nc
integer id,nts,ntc,ntt
real,dimension(:) :: Vx,Vy,Vz
real,dimension(:),optional :: T

integer numt

call nfseis_diminfo(fnm_nc,'time',numt)
if (numt<1) then
   print *, trim(fnm_nc), ' no data'
   stop 1
end if
if (numt-nts+1<ntc/ntt) then
   print *, trim(fnm_nc), ' time.length-nts+1<ntc/ntt'
   print *, numt,nts,ntc,ntt
   stop 1
end if

call nfseis_varget(fnm_nc,'Vx',Vx(1:ntc),(/id,nts/),(/1,ntc/),(/1,ntt/))
call nfseis_varget(fnm_nc,'Vy',Vy(1:ntc),(/id,nts/),(/1,ntc/),(/1,ntt/))
call nfseis_varget(fnm_nc,'Vz',Vz(1:ntc),(/id,nts/),(/1,ntc/),(/1,ntt/))
if (present(T)) then
call nfseis_varget(fnm_nc,'time',T(1:ntc),(/nts/),(/ntc/),(/ntt/))
end if
end subroutine retrieve_recvline

subroutine retrieve_snap_seis(fnm_prefix,i,j,k,varnm,var,num_nt)
character (len=*) :: fnm_prefix,varnm
integer i,j,k,num_nt
real,dimension(:) :: var
character (len=SEIS_STRLEN) :: filenm

integer m,n,mt
m=0; n=0
do 
   if (m>=num_nt) exit
   n=n+1
   filenm=trim(fnm_prefix)               &
       //'_n'//trim(io_out_pattern(n,5)) &
       //'.nc'
   call nfseis_diminfo(filenm,'time',mt)
   if (m+mt>num_nt) mt=num_nt-m
   call nfseis_varget(filenm,varnm,var(m+1:m+mt), &
        (/i,j,k,1/),(/1,1,1,mt/),(/1,1,1,1/))
   m=m+mt
end do
end subroutine retrieve_snap_seis
subroutine retrieve_snap_time(fnm_prefix,T,num_nt)
character (len=*) :: fnm_prefix
real,dimension(:) :: T
integer num_nt
character (len=SEIS_STRLEN) :: filenm

integer m,n,mt
m=0; n=0
do 
   if (m>=num_nt) exit
   n=n+1
   filenm=trim(fnm_prefix)               &
       //'_n'//trim(io_out_pattern(n,5)) &
       //'.nc'
   call nfseis_diminfo(filenm,'time',mt)
   if (m+mt>num_nt) mt=num_nt-m
   call nfseis_varget(filenm,'time',T(m+1:m+mt), &
        (/1/),(/mt/),(/1/))
   m=m+mt
end do
end subroutine retrieve_snap_time

!---------------------------------------------------------------------------

function get_fnm_seismoinfo(pnm_info,n_i,n_j,n_k) result(filenm)
  character (len=*) :: pnm_info
  integer n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_info)                         &
       //'/'//'seismoinfo'                      &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'.nc'
end function get_fnm_seismoinfo
function get_fnm_seismo(pnm,n_i,n_j,n_k) result(filenm)
  character (len=*) :: pnm
  integer n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm    )                          &
       //'/'//'seismo'                          &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'.nc'
end function get_fnm_seismo
function get_fnm_snapnode(pnm,prefix,n,ntime,n_i,n_j,n_k) result(filenm)
  integer n,ntime,n_i,n_j,n_k
  character (len=*) :: pnm,prefix
  character (len=SEIS_STRLEN) :: filenm
  integer n0
  n0=(ntime+snap_tinv(n)*snap_tcnt(n)-1)/(snap_tinv(n)*snap_tcnt(n))
  filenm=trim(pnm)                              &
       //'/'//trim(io_enum(prefix,snap_id(n)))  &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'_n'//trim(io_out_pattern(n0,5))       &
       //'.nc'
end function get_fnm_snapnode
function get_fnm_snapnode_n(pnm,prefix,id,n0,n_i,n_j,n_k) result(filenm)
  character (len=*) :: pnm,prefix
  integer id,n0,n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm)                              &
       //'/'//trim(io_enum(prefix,id))          &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k)) &
       //'_n'//trim(io_out_pattern(n0,5))       &
       //'.nc'
end function get_fnm_snapnode_n

function get_fnm_recv() result(filenm)
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_seismo) &
       //'/'//trim(fnm_recv)
end function get_fnm_recv
function get_fnm_line(id) result(filenm)
  integer id
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_seismo)                   &
       //'/'//trim(io_enum(prefix_line,id)) &
       //'.nc'
end function get_fnm_line
function get_fnm_snap(n,ntime) result(filenm)
  integer n,ntime
  character (len=SEIS_STRLEN) :: filenm
  integer n0
  n0=(ntime+snap_tinv(n)*snap_tcnt(n)-1)/(snap_tinv(n)*snap_tcnt(n))
  filenm=trim(pnm_snap)                             &
       //'/'//trim(io_enum(prefix_snap,snap_id(n))) &
       //'_n'//trim(io_out_pattern(n0,5))           &
       //'.nc'
end function get_fnm_snap
function get_fnm_snap_n(n,n0) result(filenm)
  integer n,n0
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_snap)                             &
       //'/'//trim(io_enum(prefix_snap,snap_id(n))) &
       //'_n'//trim(io_out_pattern(n0,5))           &
       //'.nc'
end function get_fnm_snap_n
function get_fnm_snapcoord(n) result(filenm)
  integer n
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_snap)                             &
       //'/'//trim(io_enum(prefix_snap,snap_id(n))) &
       //'_coord'                                   &
       //'.nc'
end function get_fnm_snapcoord

!function set_mpi_subfix(i,j,k) result(filenm)
!    integer i,j,k
!    character (len=SEIS_STRLEN) :: filenm
!    character (len=SEIS_STRLEN) :: str1,str2,str3
!    write(str1,"(i2.2)") i
!    write(str2,"(i2.2)") j
!    write(str3,"(i2.2)") k
!    filenm  ='mpi'//trim(adjustl(str1))  &
!                    //trim(adjustl(str2))  &
!                    //trim(adjustl(str3))
!end function set_mpi_subfix

function io_enum(prefix,num) result(ioname)
character (len=*) :: prefix
integer num
character (len=SEIS_STRLEN) :: ioname
character (len=SEIS_STRLEN) :: str
write(str,"(i3.3)") num
ioname=trim(prefix)//trim(str)
end function io_enum

function io_out_pattern(num,width) result(ioname)
integer num
integer,optional :: width
character (len=SEIS_STRLEN) :: ioname,str,fmt_str
if (present(width)) then
   write(str,"(i)") width
   fmt_str="(i"//trim(str)//"."//trim(str)//")"
else
   fmt_str="(i4.4)"
end if
write(ioname,fmt_str) num
end function io_out_pattern

subroutine io_delete(filenm)
character (len=*) :: filenm
integer fid
fid=5001
open(fid,file=trim(filenm),status='unknown')
close(fid,status='delete')
end subroutine io_delete

end module io_mod
