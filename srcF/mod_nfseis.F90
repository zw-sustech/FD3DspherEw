!******************************************************************************!
!*  This module relates with netcdf subroutines                               *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

#include "mod_macdrp.h"
!-----------------------------------------------------------------------------
module nfseis_mod
!-----------------------------------------------------------------------------

use netcdf
use constants_mod
implicit none

interface nfseis_attget
   module procedure nfseis_attget_int
   module procedure nfseis_attget_int1d
   module procedure nfseis_attget_real
   module procedure nfseis_attget_real1d
end interface

interface nfseis_attput
   module procedure nfseis_attput_int
   module procedure nfseis_attput_int1d
   module procedure nfseis_attput_real
   module procedure nfseis_attput_real1d
end interface

interface nfseis_varget
   module procedure nfseis_varget_int
   module procedure nfseis_varget_int1d
   module procedure nfseis_varget_int2d
   module procedure nfseis_varget_int3d
   module procedure nfseis_varget_int4d
   module procedure nfseis_varget_real
   module procedure nfseis_varget_real1d
   module procedure nfseis_varget_real2d
   module procedure nfseis_varget_real3d
   module procedure nfseis_varget_real4d
end interface

interface nfseis_varput
   module procedure nfseis_varput_int
   module procedure nfseis_varput_int1d
   module procedure nfseis_varput_int2d
   module procedure nfseis_varput_int3d
   module procedure nfseis_varput_real
   module procedure nfseis_varput_real1d
   module procedure nfseis_varput_real2d
   module procedure nfseis_varput_real3d
   module procedure nfseis_varput_real4d
end interface

interface nfseis_put
   module procedure nfseis_put_int
   module procedure nfseis_put_int1d
   module procedure nfseis_put_int2d
   module procedure nfseis_put_int3d
   module procedure nfseis_put_real
   module procedure nfseis_put_real1d
   module procedure nfseis_put_real2d
   module procedure nfseis_put_real3d
end interface

!-------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------!

!*************************************************************************
!* save seismogram on receivers                                          *
!*************************************************************************
subroutine nfseis_seismoinfo_create( filenm,tinv,title )
character (len=*) :: filenm
integer tinv
character (len=*),optional :: title
integer :: ierr,ncid,gdimid,ndimid,oldMode
integer :: indxid,gindxid,coordid,loctid,idid,twoelemid
!--
ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_handle_err(ierr)
!--
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_handle_err(ierr)
!-- define dim
ierr=nf90_def_dim(ncid, "num_pt", nf90_unlimited, ndimid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, "geo_dim", SEIS_GEO, gdimid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, "twoelem", 2, twoelemid)
     call nfseis_handle_err(ierr)
! -- define variable
ierr=nf90_def_var(ncid, 'indx', nf90_int, (/ gdimid, ndimid /), indxid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'gindx', nf90_int, (/ gdimid, ndimid /), gindxid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'coord', SEISNC_DATATYPE, (/ gdimid,ndimid /), coordid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'grid', SEISNC_DATATYPE, (/ gdimid,ndimid /), loctid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'id', nf90_int, (/ twoelemid,ndimid /), idid )
     call nfseis_handle_err(ierr)
!-- define global attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"tinv",tinv )
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
end if
!--
ierr=nf90_enddef(ncid)
!-- 
ierr=nf90_close(ncid);call nfseis_handle_err(ierr)
end subroutine nfseis_seismoinfo_create

subroutine nfseis_seismo_init(filenm,      &! file name
           num_pt,ncid,tid,vxid,vyid,vzid, &! output ids
           xid,yid,zid,                    &
           title)  ! input att
character (len=*) :: filenm
integer,intent(in) :: num_pt
integer,intent(out) :: ncid,tid,vxid,vyid,vzid
integer,intent(out),optional :: xid,yid,zid
character (len=*),intent(in),optional :: title
integer ierr,ndimid,tdimid,oldMode
!--
ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_handle_err(ierr)
!--
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_handle_err(ierr)
!-- define dim
ierr=nf90_def_dim(ncid, "num_pt", num_pt, ndimid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, "time", nf90_unlimited, tdimid)
     call nfseis_handle_err(ierr)
! -- define variable
ierr=nf90_def_var(ncid, 'time', SEISNC_DATATYPE, (/ tdimid /), tid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'Vx', SEISNC_DATATYPE, (/ ndimid, tdimid /), vxid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'Vy', SEISNC_DATATYPE, (/ ndimid, tdimid /), vyid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'Vz', SEISNC_DATATYPE, (/ ndimid, tdimid /), vzid )
     call nfseis_handle_err(ierr)
if (present(xid)) then
ierr=nf90_def_var(ncid, 'x', SEISNC_DATATYPE, (/ ndimid /), xid )
     call nfseis_handle_err(ierr)
end if
if (present(yid)) then
ierr=nf90_def_var(ncid, 'y', SEISNC_DATATYPE, (/ ndimid /), yid )
     call nfseis_handle_err(ierr)
end if
if (present(zid)) then
ierr=nf90_def_var(ncid, 'z', SEISNC_DATATYPE, (/ ndimid /), zid )
     call nfseis_handle_err(ierr)
end if
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
end if
!--
ierr=nf90_enddef(ncid)
!-- put data --
end subroutine nfseis_seismo_init
subroutine nfseis_seismo_reinit(filenm, &! file name
           ncid,tid,vxid,vyid,vzid,     &! output ids
           xid,yid,zid)  ! input att
character (len=*) :: filenm
integer,intent(out),optional :: ncid,tid,vxid,vyid,vzid
integer,intent(out),optional :: xid,yid,zid
integer ierr
!--
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
! -- inquire variable
ierr=nf90_inq_varid(ncid, 'time',  tid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, 'Vx', vxid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, 'Vy', vyid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, 'Vz', vzid )
     call nfseis_handle_err(ierr)
if (present(xid)) then
ierr=nf90_inq_varid(ncid, 'x', xid )
     call nfseis_handle_err(ierr)
end if
if (present(yid)) then
ierr=nf90_inq_varid(ncid, 'y', yid )
     call nfseis_handle_err(ierr)
end if
if (present(zid)) then
ierr=nf90_inq_varid(ncid, 'z', zid )
     call nfseis_handle_err(ierr)
end if
!--
end subroutine nfseis_seismo_reinit

!*************************************************************************
!* save wave field in 3d grid points                                     *
!*************************************************************************
subroutine nfseis_wave_init(filenm,     &
           nmvx,nmvy,nmvz,              &
           ncid,vxid,vyid,vzid,tid,nid, &
           stept,                       &
           gsubs,gsubc,gsubt,gsube,     &
           subs,subc,subt,sube,title)
character (len=*),intent(in) :: filenm,nmvx,nmvy,nmvz
integer ncid,vxid,vyid,vzid,tid,nid
real(SP),intent(in) :: stept
integer,dimension(SEIS_GEO) :: gsubs,gsubc,gsubt,gsube
integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
character (len=*),intent(in),optional :: title

integer ierr,oldMode
integer tdimid,xdimid,ydimid,zdimid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode); call nfseis_handle_err(ierr)
! -- define dim
ierr=nf90_def_dim(ncid, 'I', subc(1), xdimid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'J', subc(2), ydimid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'K', subc(3), zdimid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'time', nf90_unlimited, tdimid)
     call nfseis_handle_err(ierr)
! -- define variable
ierr=nf90_def_var(ncid,nmvx,SEISNC_DATATYPE,(/ xdimid, ydimid, zdimid,tdimid /), vxid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid,nmvy,SEISNC_DATATYPE,(/ xdimid, ydimid, zdimid,tdimid /), vyid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid,nmvz,SEISNC_DATATYPE,(/ xdimid, ydimid, zdimid,tdimid /), vzid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid,'time',SEISNC_DATATYPE,(/ tdimid /), tid )
     call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid,'ntime',nf90_int,(/ tdimid /), nid )
     call nfseis_handle_err(ierr)
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
end if
!-- define variable attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsubs",gsubs)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsubc",gsubc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsubt",gsubt)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"gsube",gsube)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"subs",subs)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"subc",subc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"subt",subt)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"sube",sube)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"stept",stept )
!--
ierr=nf90_enddef(ncid)
!--
!ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine nfseis_wave_init 
subroutine nfseis_wave_reinit(filenm,   &
           nmvx,nmvy,nmvz,              &
           ncid,vxid,vyid,vzid,tid,nid)
character (len=*),intent(in) :: filenm,nmvx,nmvy,nmvz
integer ncid,vxid,vyid,vzid,tid,nid

integer ierr

ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
! -- define variable
ierr=nf90_inq_varid(ncid, 'ntime', nid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, 'time',  tid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, nmvx, vxid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, nmvy, vyid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, nmvz, vzid )
     call nfseis_handle_err(ierr)
!--
end subroutine nfseis_wave_reinit 

!----------------------------------------------------------------

subroutine nfseis_data_create(filenm, &
                       nx,ny,nz,      &
                       title)
character (len=*),intent(in) :: filenm
integer :: nx,ny,nz
character (len=*),intent(in),optional :: title

integer ncid,ierr,oldMode
integer nxid,nyid,nzid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode);    call nfseis_handle_err(ierr)
! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, nxid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'J', ny, nyid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'K', nz, nzid);     call nfseis_handle_err(ierr)
!-- define global attribute
if (present(title) ) ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
!--
ierr=nf90_enddef(ncid); call nfseis_handle_err(ierr)
!--
ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine nfseis_data_create

subroutine nfseis_data_addvar(filenm,aVarNm)
character (len=*),intent(in) :: aVarNm,filenm

integer ncid,ierr,oldMode
integer xdimid,ydimid,zdimid
integer varid

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid );call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid);                          call nfseis_handle_err(ierr)
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode);    call nfseis_handle_err(ierr)
! -- req dim
ierr=nf90_inq_dimid(ncid,'I',xdimid); call nfseis_handle_err(ierr)
ierr=nf90_inq_dimid(ncid,'J',ydimid); call nfseis_handle_err(ierr)
ierr=nf90_inq_dimid(ncid,'K',zdimid); call nfseis_handle_err(ierr)
! -- define variable
ierr=nf90_def_var(ncid, trim(aVarNm), SEISNC_DATATYPE, &
                 (/ xdimid, ydimid, zdimid /),    &
                 varid )
                 call nfseis_handle_err(ierr)
!--
ierr=nf90_enddef(ncid)
!--
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_data_addvar

subroutine nfseis_data_attput(filenm,         &
           subs,subc,subt,                    &
           ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,  &
           nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz,  &
           ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,     &
           ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,     &
           point_in_this)
character (len=*),intent(in) :: filenm
integer,dimension(SEIS_GEO) :: subs,subc,subt
integer,intent(in) ::                         &
        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,     &
        nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz,     &
        ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,        &
        ngx1,ngx2,ngy1,ngy2,ngz1,ngz2 !,      &
integer,dimension(SEIS_GEO*2) :: point_in_this

integer ierr,ncid

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid );call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid);                          call nfseis_handle_err(ierr)
!-- define variable attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ni1",ni1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ni2",ni2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nj1",nj1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nj2",nj2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nk1",nk1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nk2",nk2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ni",ni)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nj",nj)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nk",nk)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nx1",nx1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nx2",nx2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ny1",ny1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ny2",ny2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nz1",nz1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nz2",nz2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nx",nx)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ny",ny)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"nz",nz)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngi1",ngi1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngi2",ngi2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngj1",ngj1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngj2",ngj2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngk1",ngk1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngk2",ngk2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngx1",ngx1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngx2",ngx2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngy1",ngy1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngy2",ngy2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngz1",ngz1)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"ngz2",ngz2)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"point_in_this",point_in_this)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"start",subs)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"count",subc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"stride",subt)
!--
ierr=nf90_enddef(ncid); call nfseis_handle_err(ierr)
!--
ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine nfseis_data_attput

subroutine nfseis_data_attget(filenm,        &
           ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &
           nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &
           ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &
           ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,    &
           point_in_this                     &
           )
character (len=*),intent(in) :: filenm
integer,dimension(SEIS_GEO*2),optional :: point_in_this
integer,optional,intent(out) ::           &
        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &
        nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &
        ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &
        ngx1,ngx2,ngy1,ngy2,ngz1,ngz2 !, &

integer ncid,ierr

ierr=nf90_open(trim(filenm),IOR(NF90_NOWRITE,NF90_SHARE),ncid )
call nfseis_handle_err(ierr)
!IOR(NF_WRITE, NF_SHARE)
if (present(nx1))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nx1',nx1 )
if (present(nx2))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nx2',nx2 )
if (present(ny1))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'ny1',ny1 )
if (present(ny2))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'ny2',ny2 )
if (present(nz1))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nz1',nz1 )
if (present(nz2))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nz2',nz2 )
if (present(nx))   ierr=nf90_get_att(ncid,NF90_GLOBAL,'nx', nx )
if (present(ny))   ierr=nf90_get_att(ncid,NF90_GLOBAL,'ny', ny )
if (present(nz))   ierr=nf90_get_att(ncid,NF90_GLOBAL,'nz', nz )
if (present(ni1))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'ni1',ni1 )
if (present(ni2))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'ni2',ni2 )
if (present(nj1))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nj1',nj1 )
if (present(nj2))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nj2',nj2 )
if (present(nk1))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nk1',nk1 )
if (present(nk2))  ierr=nf90_get_att(ncid,NF90_GLOBAL,'nk2',nk2 )
if (present(ni))   ierr=nf90_get_att(ncid,NF90_GLOBAL,'ni', ni )
if (present(nj))   ierr=nf90_get_att(ncid,NF90_GLOBAL,'nj', nj )
if (present(nk))   ierr=nf90_get_att(ncid,NF90_GLOBAL,'nk', nk )
if (present(ngi1)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngi1', ngi1 )
if (present(ngi2)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngi2', ngi2 )
if (present(ngj1)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngj1', ngj1 )
if (present(ngj2)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngj2', ngj2 )
if (present(ngk1)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngk1', ngk1 )
if (present(ngk2)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngk2', ngk2 )
if (present(ngx1)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngx1', ngx1 )
if (present(ngx2)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngx2', ngx2 )
if (present(ngy1)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngy1', ngy1 )
if (present(ngy2)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngy2', ngy2 )
if (present(ngz1)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngz1', ngz1 )
if (present(ngz2)) ierr=nf90_get_att(ncid,NF90_GLOBAL,'ngz2', ngz2 )
if (present(point_in_this)) &
    ierr=nf90_get_att(ncid,NF90_GLOBAL,'point_in_this',point_in_this )
ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine nfseis_data_attget

!----------------------------------------------------------------

subroutine nfseis_srcnode_create(filenm,   &
        ntwin,nslen)
character (len=*) :: filenm
integer ntwin,nslen

integer ncid,ierr,oldMode
integer npid,ntid,twoid,geoid,slenid
integer txxid,tyyid,tzzid,txyid,txzid,tyzid
integer RuptTid,indxid,stfaid,twinid,sdeltid,siftid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode); call nfseis_handle_err(ierr)

! -- define dim
ierr=nf90_def_dim(ncid,'npt',nf90_unlimited ,npid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid,'ntwin' , ntwin ,ntid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid,'Twoelem',2,twoid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid,'geodim',SEIS_GEO,geoid)
     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid,'nslen',nslen,slenid)
     call nfseis_handle_err(ierr)

! -- define variable
ierr=nf90_def_var(ncid, 'MomTxx', SEISNC_DATATYPE,     &
                  (/ ntid, npid /),               &
                  txxid )
                  call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'MomTyy', SEISNC_DATATYPE,     &
                  (/ ntid, npid /),               &
                 tyyid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'MomTzz', SEISNC_DATATYPE,     &
                  (/ ntid, npid /),               &
                 tzzid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'MomTxy', SEISNC_DATATYPE,     &
                  (/ ntid, npid /),               &
                 txyid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'MomTxz', SEISNC_DATATYPE,     &
                  (/ ntid, npid /),               &
                 txzid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'MomTyz', SEISNC_DATATYPE,     &
                  (/ ntid, npid /),               &
                 tyzid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'RuptT', SEISNC_DATATYPE,      &
                  (/ npid /),                     &
                 RuptTid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'fault_indx', nf90_int,   &
                  (/ geoid, npid /),              &
                 indxid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'fault_loct', nf90_float,   &
                  (/ geoid, npid /),              &
                 indxid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'sdelt', SEISNC_DATATYPE,      &
                 (/ slenid,slenid,slenid,npid /), &
                 sdeltid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'shift', SEISNC_DATATYPE,      &
                 (/ geoid,npid /), &
                 siftid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'stf_A', SEISNC_DATATYPE,      &
                  (/ ntid /),                     &
                 stfaid )
                 call nfseis_handle_err(ierr)
ierr=nf90_def_var(ncid, 'twin', SEISNC_DATATYPE,       &
                  (/ twoid, ntid /),              &
                 twinid )
                 call nfseis_handle_err(ierr)
ierr=nf90_enddef(ncid)
!-- close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_srcnode_create

subroutine nfseis_srcnode_attput(filenm,   &
        flag_stf_type, flag_mech_type, stf_t0, stf_alpha)
character (len=*) :: filenm
real(SP) :: stf_t0,stf_alpha
integer :: flag_stf_type,flag_mech_type

integer ncid,ierr

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid)
     call nfseis_handle_err(ierr)

!-- define global attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"flag_stf_type",flag_stf_type)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"flag_mech_type",flag_mech_type)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"stf_t0",stf_t0)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"stf_alpha",stf_alpha)

!--
ierr=nf90_enddef(ncid)
!-- close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_srcnode_attput

subroutine nfseis_grid_create_base(filenm,nx,ny,nz)
character (len=*),intent(in) :: filenm
integer :: nx,ny,nz

integer ncid,ierr,oldMode
integer xid,yid,zid
integer coordxid,coordyid,coordzid,zhalfid,zgridid

ierr=nf90_create( path=trim(filenm),cmode=nf90_clobber,ncid=ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_handle_err(ierr)

! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, xid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'J', ny, yid);     call nfseis_handle_err(ierr)
ierr=nf90_def_dim(ncid, 'K', nz, zid);     call nfseis_handle_err(ierr)

! -- variable
ierr=nf90_def_var(ncid, 'x', nf90_float, (/ xid /), coordxid )
ierr=nf90_def_var(ncid, 'y', nf90_float, (/ yid /), coordyid )
ierr=nf90_def_var(ncid, 'z', nf90_float, (/ zid /), coordzid )
ierr=nf90_def_var(ncid, 'xsin', nf90_float, (/ xid /), coordyid )
ierr=nf90_def_var(ncid, 'xcot', nf90_float, (/ xid /), coordyid )
!--
ierr=nf90_enddef(ncid)
!--
ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine nfseis_grid_create_base

!----------------------------------------------------------------

subroutine nfseis_diminfo(filenm,dimnm,dimlen)
character (len=*) :: filenm,dimnm
integer,optional :: dimlen

integer :: ierr,ncid,dimid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
! get diminfo
ierr=nf90_inq_dimid(ncid, dimnm, dimid )
     call nfseis_handle_err(ierr)
ierr=nf90_inquire_dimension( ncid, dimid, len=dimlen)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_diminfo

! --- put ---
subroutine nfseis_put_int(ncid,varid,var,subs,subc,subt)
integer ncid,varid
integer :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_int
subroutine nfseis_put_int1d(ncid,varid,var,subs,subc,subt)
integer ncid,varid
integer,dimension(:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_int1d
subroutine nfseis_put_int2d(ncid,varid,var,subs,subc,subt)
integer ncid,varid
integer,dimension(:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_int2d
subroutine nfseis_put_int3d(ncid,varid,var,subs,subc,subt)
integer ncid,varid
integer,dimension(:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_int3d
subroutine nfseis_put_real(ncid,varid,var,subs,subc,subt)
integer ncid,varid
real(SP) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_real
subroutine nfseis_put_real1d(ncid,varid,var,subs,subc,subt)
integer ncid,varid
real(SP),dimension(:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_real1d
subroutine nfseis_put_real2d(ncid,varid,var,subs,subc,subt)
integer ncid,varid
real(SP),dimension(:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_real2d
subroutine nfseis_put_real3d(ncid,varid,var,subs,subc,subt)
integer ncid,varid
real(SP),dimension(:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
end subroutine nfseis_put_real3d

! --- varput ---
subroutine nfseis_varput_int(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_int
subroutine nfseis_varput_int1d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_int1d
subroutine nfseis_varput_int2d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_int2d
subroutine nfseis_varput_int3d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_int3d
subroutine nfseis_varput_real(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_real
subroutine nfseis_varput_real1d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_real1d
subroutine nfseis_varput_real2d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_real2d
subroutine nfseis_varput_real3d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_real3d
subroutine nfseis_varput_real4d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:,:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varput_real4d

! --- get var ---
subroutine nfseis_varget_int(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
integer,dimension(1:1) :: a
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,a,subs,subc,subt)
     call nfseis_handle_err(ierr)
var=a(1)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_int
subroutine nfseis_varget_int1d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_int1d
subroutine nfseis_varget_int2d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_int2d
subroutine nfseis_varget_int3d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_int3d
subroutine nfseis_varget_int4d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
integer,dimension(:,:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_int4d
subroutine nfseis_varget_real(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
real(SP),dimension(1:1) :: a
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,a,subs,subc,subt)
     call nfseis_handle_err(ierr)
var=a(1)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_real
subroutine nfseis_varget_real1d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_real1d
subroutine nfseis_varget_real2d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_real2d
subroutine nfseis_varget_real3d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_real3d
subroutine nfseis_varget_real4d(filenm,varnm,var,subs,subc,subt)
character (len=*) :: filenm,varnm
real(SP),dimension(:,:,:,:) :: var
integer,dimension(:) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_inq_varid(ncid, varnm, varid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_varget_real4d

! --- get att ---
subroutine nfseis_attget_int(filenm,attnm,att)
character (len=*) :: filenm,attnm
integer :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
   call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attget_int
subroutine nfseis_attget_int1d(filenm,attnm,att)
character (len=*) :: filenm,attnm
integer,dimension(:) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
   call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attget_int1d
subroutine nfseis_attget_real(filenm,attnm,att)
character (len=*) :: filenm,attnm
real(SP) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
   call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attget_real
subroutine nfseis_attget_real1d(filenm,attnm,att)
character (len=*) :: filenm,attnm
real(SP),dimension(:) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
   call nfseis_handle_err(ierr)
! close
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attget_real1d

! --- put att ---
subroutine nfseis_attput_int(filenm,attnm,att)
character (len=*) :: filenm,attnm
integer :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid)
     call nfseis_handle_err(ierr)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_enddef(ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attput_int
subroutine nfseis_attput_int1d(filenm,attnm,att)
character (len=*) :: filenm,attnm
integer,dimension(:) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid)
     call nfseis_handle_err(ierr)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_enddef(ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attput_int1d
subroutine nfseis_attput_real(filenm,attnm,att)
character (len=*) :: filenm,attnm
real(SP) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid)
     call nfseis_handle_err(ierr)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_enddef(ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attput_real
subroutine nfseis_attput_real1d(filenm,attnm,att)
character (len=*) :: filenm,attnm
real(SP),dimension(:) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_redef( ncid)
     call nfseis_handle_err(ierr)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
     call nfseis_handle_err(ierr)
! close
ierr=nf90_enddef(ncid)
     call nfseis_handle_err(ierr)
ierr=nf90_close(ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_attput_real1d

subroutine nfseis_open(filenm,ncid)
character (len=*) :: filenm
integer :: ierr,ncid
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_handle_err(ierr)
end subroutine nfseis_open
subroutine nfseis_inq_varid(ncid,vnm,vid)
integer,intent(in) :: ncid
integer,intent(out) :: vid
character (len=*) :: vnm
integer ierr
ierr=nf90_inq_varid(ncid, vnm, vid )
     call nfseis_handle_err(ierr)
end subroutine nfseis_inq_varid
subroutine nfseis_close(ncid)
    integer ncid,ierr
    ierr=nf90_close(ncid);call nfseis_handle_err(ierr)
end subroutine nfseis_close

subroutine nfseis_handle_err(ierr)
    integer,intent(in) :: ierr
    if (ierr /= nf90_noerr) then
       print *, trim(nf90_strerror(ierr))
       stop 2
    end if
end subroutine nfseis_handle_err

end module nfseis_mod

