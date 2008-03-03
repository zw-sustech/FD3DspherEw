!******************************************************************************!
!*  This module relates with netcdf subroutines                               *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
module nfseis_mod
!-----------------------------------------------------------------------------

#include "mod_macdrp.h"

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
character (len=*),intent(in) :: filenm
integer,intent(in) :: tinv
character (len=*),optional,intent(in) :: title
integer :: ierr,ncid,gdimid,ndimid,oldMode
integer :: indxid,gindxid,coordid,loctid,idid,twoelemid
!--
ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'seismoinfo_create:'//trim(filenm))
!--
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in seismoinfo_create')
!-- define dim
ierr=nf90_def_dim(ncid, "num_pt", nf90_unlimited, ndimid)
     call nfseis_except(ierr,'num_pt dim in seismoinfo_create')
ierr=nf90_def_dim(ncid, "geo_dim", SEIS_GEO, gdimid)
     call nfseis_except(ierr,'geo_dim dim in seismoinfo_create')
ierr=nf90_def_dim(ncid, "twoelem", 2, twoelemid)
     call nfseis_except(ierr,'twoelem dim in seismoinfo_create')
! -- define variable
ierr=nf90_def_var(ncid, 'indx', nf90_int, (/ gdimid, ndimid /), indxid )
     call nfseis_except(ierr,'indx var in seismoinfo_create')
ierr=nf90_def_var(ncid, 'gindx', nf90_int, (/ gdimid, ndimid /), gindxid )
     call nfseis_except(ierr,'gindx var in seismoinfo_create')
ierr=nf90_def_var(ncid, 'coord', SEISNC_DATATYPE, (/ gdimid,ndimid /), coordid )
     call nfseis_except(ierr,'coord var in seismoinfo_create')
ierr=nf90_def_var(ncid, 'grid', SEISNC_DATATYPE, (/ gdimid,ndimid /), loctid )
     call nfseis_except(ierr,'grid var in seismoinfo_create')
ierr=nf90_def_var(ncid, 'id', nf90_int, (/ twoelemid,ndimid /), idid )
     call nfseis_except(ierr,'id var in seismoinfo_create')
!-- define global attribute
ierr=nf90_put_att(ncid,NF90_GLOBAL,"tinv",tinv )
     call nfseis_except(ierr,'tinv att in seismoinfo_create')
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in seismoinfo_create')
end if
!--
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in seismoinfo_create')
!-- 
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in seismoinfo_create')
end subroutine nfseis_seismoinfo_create

subroutine nfseis_seismo_init(filenm,      &
           num_pt,ncid,tid,vxid,vyid,vzid, &
           xid,yid,zid,                    &
           title)  ! input att
character (len=*),intent(in) :: filenm
integer,intent(in) :: num_pt
integer,intent(out) :: ncid,tid,vxid,vyid,vzid
integer,intent(out),optional :: xid,yid,zid
character (len=*),intent(in),optional :: title
integer ierr,ndimid,tdimid,oldMode
!--
ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'seismo_init:'//trim(filenm))
!--
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in seismo_init')
!-- define dim
ierr=nf90_def_dim(ncid, "num_pt", num_pt, ndimid)
     call nfseis_except(ierr,'num_pt dim in seismo_init')
ierr=nf90_def_dim(ncid, "time", nf90_unlimited, tdimid)
     call nfseis_except(ierr,'time dim in seismo_init')
! -- define variable
ierr=nf90_def_var(ncid, 'time', SEISNC_DATATYPE, (/ tdimid /), tid )
     call nfseis_except(ierr,'time var in seismo_init')
ierr=nf90_def_var(ncid, 'Vx', SEISNC_DATATYPE, (/ ndimid, tdimid /), vxid )
     call nfseis_except(ierr,'Vx var in seismo_init')
ierr=nf90_def_var(ncid, 'Vy', SEISNC_DATATYPE, (/ ndimid, tdimid /), vyid )
     call nfseis_except(ierr,'Vy var in seismo_init')
ierr=nf90_def_var(ncid, 'Vz', SEISNC_DATATYPE, (/ ndimid, tdimid /), vzid )
     call nfseis_except(ierr,'Vz var in seismo_init')
if (present(xid)) then
ierr=nf90_def_var(ncid, 'x', SEISNC_DATATYPE, (/ ndimid /), xid )
     call nfseis_except(ierr,'x var in seismo_init')
end if
if (present(yid)) then
ierr=nf90_def_var(ncid, 'y', SEISNC_DATATYPE, (/ ndimid /), yid )
     call nfseis_except(ierr,'y var in seismo_init')
end if
if (present(zid)) then
ierr=nf90_def_var(ncid, 'z', SEISNC_DATATYPE, (/ ndimid /), zid )
     call nfseis_except(ierr,'z var in seismo_init')
end if
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in seismo_init')
end if
!--
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in seismo_init')
!-- put data --
end subroutine nfseis_seismo_init
subroutine nfseis_seismo_reinit(filenm, &! file name
           ncid,tid,vxid,vyid,vzid,     &! output ids
           xid,yid,zid)  ! input att
character (len=*),intent(in) :: filenm
integer,intent(out),optional :: ncid,tid,vxid,vyid,vzid
integer,intent(out),optional :: xid,yid,zid
integer ierr
!--
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     call nfseis_except(ierr,'seismo_reinit:'//trim(filenm))
! -- inquire variable
ierr=nf90_inq_varid(ncid, 'time',  tid )
     call nfseis_except(ierr,'time var in seismo_reinit')
ierr=nf90_inq_varid(ncid, 'Vx', vxid )
     call nfseis_except(ierr,'Vx var in seismo_reinit')
ierr=nf90_inq_varid(ncid, 'Vy', vyid )
     call nfseis_except(ierr,'Vy var in seismo_reinit')
ierr=nf90_inq_varid(ncid, 'Vz', vzid )
     call nfseis_except(ierr,'Vz var in seismo_reinit')
if (present(xid)) then
ierr=nf90_inq_varid(ncid, 'x', xid )
     call nfseis_except(ierr,'x var in seismo_reinit')
end if
if (present(yid)) then
ierr=nf90_inq_varid(ncid, 'y', yid )
     call nfseis_except(ierr,'y var in seismo_reinit')
end if
if (present(zid)) then
ierr=nf90_inq_varid(ncid, 'z', zid )
     call nfseis_except(ierr,'z var in seismo_reinit')
end if
!--
end subroutine nfseis_seismo_reinit

!*************************************************************************
!* save wave field in 3d grid points                                     *
!*************************************************************************
subroutine nfseis_snap_header(filenm,     &
           ncid,tid, &
           stept,                       &
           gsubs,gsubc,gsubt,gsube,     &
           subs,subc,subt,sube,title)
character (len=*),intent(in) :: filenm
integer,intent(out) :: ncid,tid
real(SP),intent(in) :: stept
integer,dimension(SEIS_GEO),intent(in) :: gsubs,gsubc,gsubt,gsube
integer,dimension(SEIS_GEO),intent(in) :: subs,subc,subt,sube
character (len=*),intent(in),optional :: title

integer ierr,oldMode
integer tdimid,xdimid,ydimid,zdimid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'wave_init:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in wave_init')
! -- define dim
ierr=nf90_def_dim(ncid, 'I', subc(1), xdimid)
     call nfseis_except(ierr,'I dim in wave_init')
ierr=nf90_def_dim(ncid, 'J', subc(2), ydimid)
     call nfseis_except(ierr,'J dim in wave_init')
ierr=nf90_def_dim(ncid, 'K', subc(3), zdimid)
     call nfseis_except(ierr,'K dim in wave_init')
ierr=nf90_def_dim(ncid, 'time', nf90_unlimited, tdimid)
     call nfseis_except(ierr,'time dim in wave_init')
ierr=nf90_def_var(ncid,'time',SEISNC_DATATYPE,(/ tdimid /), tid )
     call nfseis_except(ierr,'time var in wave_init')
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in wave_init')
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
!ierr=nf90_enddef(ncid)
!     call nfseis_except(ierr,'enddef in wave_init')
end subroutine nfseis_snap_header
subroutine nfseis_snap_addvar(ncid,varnm,vid)
character (len=*),intent(in) :: varnm
integer,intent(out) :: ncid,vid
integer ierr,tdimid,xdimid,ydimid,zdimid
ierr=nf90_inq_dimid(ncid,'time',tdimid)
     call nfseis_except(ierr,'time dim in snap_addvar')
ierr=nf90_inq_dimid(ncid,'I',xdimid)
     call nfseis_except(ierr,'I dim in snap_addvar')
ierr=nf90_inq_dimid(ncid,'J',ydimid)
     call nfseis_except(ierr,'J dim in snap_addvar')
ierr=nf90_inq_dimid(ncid,'K',zdimid)
     call nfseis_except(ierr,'K dim in snap_addvar')
ierr=nf90_def_var(ncid,varnm,SEISNC_DATATYPE, &
     (/ xdimid, ydimid, zdimid,tdimid /), vid )
     call nfseis_except(ierr,'add var in snap_addvar')
end subroutine nfseis_snap_addvar
subroutine nfseis_snap_enddef(ncid)
integer ncid
integer ierr
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in snap_enddef')
end subroutine nfseis_snap_enddef

!----------------------------------------------------------------

subroutine nfseis_data_create(filenm, &
                       nx,ny,nz,      &
                       title)
character (len=*),intent(in) :: filenm
integer,intent(in) :: nx,ny,nz
character (len=*),intent(in),optional :: title

integer ncid,ierr,oldMode
integer nxid,nyid,nzid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'data_create:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in data_create')
! -- define dim
ierr=nf90_def_dim(ncid, 'I', nx, nxid)
     call nfseis_except(ierr,'I dim in data_create')
ierr=nf90_def_dim(ncid, 'J', ny, nyid)
     call nfseis_except(ierr,'J dim in data_create')
ierr=nf90_def_dim(ncid, 'K', nz, nzid)
     call nfseis_except(ierr,'K dim in data_create')
!-- define global attribute
if (present(title) ) then
   ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
     call nfseis_except(ierr,'title att in data_create')
end if
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'endedf in data_create')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in data_create')
end subroutine nfseis_data_create

subroutine nfseis_data_addvar(filenm,aVarNm)
character (len=*),intent(in) :: aVarNm,filenm

integer ncid,ierr,oldMode
integer xdimid,ydimid,zdimid
integer varid

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     call nfseis_except(ierr,'data_addvar:'//trim(filenm))
ierr=nf90_redef( ncid)
     call nfseis_except(ierr,'redef in data_addvar')
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in data_addvar')
! -- req dim
ierr=nf90_inq_dimid(ncid,'I',xdimid)
     call nfseis_except(ierr,'I dim in data_addvar')
ierr=nf90_inq_dimid(ncid,'J',ydimid)
     call nfseis_except(ierr,'J dim in data_addvar')
ierr=nf90_inq_dimid(ncid,'K',zdimid)
     call nfseis_except(ierr,'K dim in data_addvar')
! -- define variable
ierr=nf90_def_var(ncid, trim(aVarNm), SEISNC_DATATYPE, &
                 (/ xdimid, ydimid, zdimid /),         &
                 varid )
     call nfseis_except(ierr,'var def in data_addvar')
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in data_addvar')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in data_addvar')
end subroutine nfseis_data_addvar

subroutine nfseis_data_attput(filenm,         &
           subs,subc,subt,                    &
           ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,  &
           nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz,  &
           ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,     &
           ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,     &
           point_in_this)
character (len=*),intent(in) :: filenm
integer,dimension(SEIS_GEO),intent(in) :: subs,subc,subt
integer,intent(in) ::                         &
        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,     &
        nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz,     &
        ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,        &
        ngx1,ngx2,ngy1,ngy2,ngz1,ngz2 !,      &
integer,dimension(SEIS_GEO*2) :: point_in_this

integer ierr,ncid

ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
ierr=nf90_redef( ncid)
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
ierr=nf90_enddef(ncid)
!--
ierr=nf90_close(ncid)
end subroutine nfseis_data_attput

!----------------------------------------------------------------

subroutine nfseis_diminfo(filenm,dimnm,dimlen)
character (len=*),intent(in) :: filenm,dimnm
integer,intent(out) :: dimlen

integer :: ierr,ncid,dimid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'diminfo:'//trim(filenm))
ierr=nf90_inq_dimid(ncid, dimnm, dimid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'dim name in diminfo:'//trim(dimnm))
ierr=nf90_inquire_dimension( ncid, dimid, len=dimlen)
     call nfseis_except(ierr,'inquire dimlen in diminfo')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in diminfo')
end subroutine nfseis_diminfo

! --- put ---
subroutine nfseis_put_int(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int')
     end if
end subroutine nfseis_put_int
subroutine nfseis_put_int1d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int1d')
     end if
end subroutine nfseis_put_int1d
subroutine nfseis_put_int2d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int2d')
     end if
end subroutine nfseis_put_int2d
subroutine nfseis_put_int3d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
integer,dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_int3d')
     end if
end subroutine nfseis_put_int3d
subroutine nfseis_put_real(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real')
     end if
end subroutine nfseis_put_real
subroutine nfseis_put_real1d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real1d')
     end if
end subroutine nfseis_put_real1d
subroutine nfseis_put_real2d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real2d')
     end if
end subroutine nfseis_put_real2d
subroutine nfseis_put_real3d(ncid,varid,var,subs,subc,subt)
integer,intent(in) :: ncid,varid
real(SP),dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'in nfseis_put_real3d')
     end if
end subroutine nfseis_put_real3d

! --- varput ---
subroutine nfseis_varput_int(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int')
end subroutine nfseis_varput_int
subroutine nfseis_varput_int1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int1d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int1d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int1d')
end subroutine nfseis_varput_int1d
subroutine nfseis_varput_int2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int2d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int2d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int2d')
end subroutine nfseis_varput_int2d
subroutine nfseis_varput_int3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_int3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_int3d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_int3d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_int3d')
end subroutine nfseis_varput_int3d
subroutine nfseis_varput_real(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,(/var/),subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real')
end subroutine nfseis_varput_real
subroutine nfseis_varput_real1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real1d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real1d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real1d')
end subroutine nfseis_varput_real1d
subroutine nfseis_varput_real2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real2d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real2d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real2d')
end subroutine nfseis_varput_real2d
subroutine nfseis_varput_real3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real3d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real3d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real3d')
end subroutine nfseis_varput_real3d
subroutine nfseis_varput_real4d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:,:),intent(in) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open( trim(filenm), NF90_WRITE, ncid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'varput_real4d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in varput_real4d:'//trim(varnm))
ierr=nf90_put_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var put in varput_real4d')
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in varput_real4d')
end subroutine nfseis_varput_real4d

! --- get var ---
subroutine nfseis_varget_int(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
integer,dimension(1:1) :: a
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,a,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int error:'//trim(filenm))
     end if
var=a(1)
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int')
end subroutine nfseis_varget_int
subroutine nfseis_varget_int1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int1d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int1d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int1d')
end subroutine nfseis_varget_int1d
subroutine nfseis_varget_int2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int2d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int2d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int2d')
end subroutine nfseis_varget_int2d
subroutine nfseis_varget_int3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int3d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int3d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int3d')
end subroutine nfseis_varget_int3d
subroutine nfseis_varget_int4d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
integer,dimension(:,:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_int4d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_int4d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_int4d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_int4d')
end subroutine nfseis_varget_int4d
subroutine nfseis_varget_real(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
real(SP),dimension(1:1) :: a
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,a,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real error:'//trim(filenm))
     end if
var=a(1)
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real')
end subroutine nfseis_varget_real
subroutine nfseis_varget_real1d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real1d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real1d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real1d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real1d')
end subroutine nfseis_varget_real1d
subroutine nfseis_varget_real2d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real2d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real2d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real2d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real2d')
end subroutine nfseis_varget_real2d
subroutine nfseis_varget_real3d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real3d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real3d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real3d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real3d')
end subroutine nfseis_varget_real3d
subroutine nfseis_varget_real4d(filenm,varnm,var,subs,subc,subt)
character (len=*),intent(in) :: filenm,varnm
real(SP),dimension(:,:,:,:),intent(out) :: var
integer,dimension(:),intent(in) :: subs,subc,subt
integer :: ierr,ncid,varid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'vargut_real4d:'//trim(filenm))
ierr=nf90_inq_varid(ncid, varnm, varid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'var name in vargut_real4d:'//trim(varnm))
ierr=nf90_get_var(ncid,varid,var,subs,subc,subt)
     if (ierr /= nf90_noerr) then
     print *, 'subs=',subs
     print *, 'subc=',subc
     print *, 'subt=',subt
     call nfseis_except(ierr,'vargut_real4d error:'//trim(filenm))
     end if
ierr=nf90_close(ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'file close in vargut_real4d')
end subroutine nfseis_varget_real4d

! --- get att ---
subroutine nfseis_attget_int(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_int
subroutine nfseis_attget_int1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,dimension(:),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_int1d
subroutine nfseis_attget_real(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_real
subroutine nfseis_attget_real1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),dimension(:),intent(out) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_NOWRITE,ncid)
ierr=nf90_get_att(ncid,NF90_GLOBAL,attnm,att )
! close
ierr=nf90_close(ncid)
end subroutine nfseis_attget_real1d

! --- put att ---
subroutine nfseis_attput_int(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_int
subroutine nfseis_attput_int1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
integer,dimension(:),intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_int1d
subroutine nfseis_attput_real(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_real
subroutine nfseis_attput_real1d(filenm,attnm,att)
character (len=*),intent(in) :: filenm,attnm
real(SP),dimension(:),intent(in) :: att
integer :: ierr,ncid
! open
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
ierr=nf90_redef( ncid)
! put att
ierr=nf90_put_att(ncid,NF90_GLOBAL,attnm,att)
! close
ierr=nf90_enddef(ncid)
ierr=nf90_close(ncid)
end subroutine nfseis_attput_real1d

subroutine nfseis_open(filenm,ncid)
character (len=*),intent(in) :: filenm
integer,intent(out) :: ncid
integer :: ierr
ierr=nf90_open(trim(filenm),NF90_WRITE,ncid)
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'nfseis_open:'//trim(filenm))
end subroutine nfseis_open
subroutine nfseis_inq_varid(ncid,vnm,vid)
integer,intent(in) :: ncid
integer,intent(out) :: vid
character (len=*),intent(in) :: vnm
integer ierr
ierr=nf90_inq_varid(ncid, vnm, vid )
     if (ierr /= nf90_noerr) &
     call nfseis_except(ierr,'nfseis_inq_varid:'//trim(vnm))
end subroutine nfseis_inq_varid
subroutine nfseis_close(ncid)
    integer,intent(in) :: ncid
    integer ierr
    ierr=nf90_close(ncid)
end subroutine nfseis_close

subroutine nfseis_except(ierr,msg)
    integer,intent(in) :: ierr
    character (len=*),intent(in) :: msg
    if (ierr /= nf90_noerr) then
       print *, trim(msg)
       print *, trim(nf90_strerror(ierr))
       stop 2
    end if
end subroutine nfseis_except
subroutine nfseis_handle_err(ierr)
    integer,intent(in) :: ierr
    if (ierr /= nf90_noerr) then
       print *, trim(nf90_strerror(ierr))
       stop 2
    end if
end subroutine nfseis_handle_err

end module nfseis_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
