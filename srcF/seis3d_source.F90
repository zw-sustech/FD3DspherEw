!******************************************************************************!
!*  This program locates the seismic source                                   *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
program seis3d_source
!-----------------------------------------------------------------------------

#include "mod_macdrp.h"

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

real(SP) :: src_hyper_height
character (len=SEIS_STRLEN) :: filenm
real,allocatable :: gx(:),gy(:),gz(:)
integer :: nloc(1)
integer :: n_i,n_j,n_k,n1,n2,n3
integer :: n,m,nfrc,nmom,gi,gj,gk
integer :: ncid,indxid,axisid,t0id,stftid,stffid
integer :: fxid,fyid,fzid,mxxid,myyid,mzzid,mxyid,mxzid,myzid

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
!force
do n=1,num_force
   nloc=minloc(abs(force_axis(1,n)-gx)); force_indx(1,n)=loct_i(nloc(1))
   nloc=minloc(abs(force_axis(2,n)-gy)); force_indx(2,n)=loct_j(nloc(1))
   nloc=minloc(abs(force_axis(3,n)-gz)); force_indx(3,n)=loct_k(nloc(1))
   if (force_axis(3,n)>src_hyper_height) force_indx(3,n)=n3-LenFD
   gi=force_indx(1,n);gj=force_indx(2,n); gk=force_indx(3,n)
   if (     gi<ni1 .or. gi>n1-LenFD  &
       .or. gj<nj1 .or. gj>n2-LenFD  &
       .or. gk<nk1 .or. gk>n3-LenFD ) then
      print *, 'force source out the model'
      print *, 'source point is ',force_axis(:,n)
      print *, 'x scale of model:',gx(ni1),gx(n1-LenFD)
      print *, 'y scale of model:',gy(nj1),gy(n2-LenFD)
      print *, 'z scale of model:',gz(nk1),gz(n3-LenFD)
      stop 1
   end if
end do
!moment
do n=1,num_moment
   nloc=minloc(abs(moment_axis(1,n)-gx)); moment_indx(1,n)=loct_i(nloc(1))
   nloc=minloc(abs(moment_axis(2,n)-gy)); moment_indx(2,n)=loct_j(nloc(1))
   nloc=minloc(abs(moment_axis(3,n)-gz)); moment_indx(3,n)=loct_k(nloc(1))
   if (moment_axis(3,n)>src_hyper_height) moment_indx(3,n)=n3-LenFD
   gi=moment_indx(1,n);gj=moment_indx(2,n); gk=moment_indx(3,n)
   if (     gi<ni1 .or. gi>n1-LenFD  &
       .or. gj<nj1 .or. gj>n2-LenFD  &
       .or. gk<nk1 .or. gk>n3-LenFD ) then
      print *, 'moment source out the model'
      print *, 'source point is ',moment_axis(:,n)
      print *, 'x scale of model:',gx(ni1),gx(n1-LenFD)
      print *, 'y scale of model:',gy(nj1),gy(n2-LenFD)
      print *, 'z scale of model:',gz(nk1),gz(n3-LenFD)
      stop 1
   end if
end do

! distribute to thread
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)

   nfrc=0
   do n=1,num_force
      gi=force_indx(1,n);gj=force_indx(2,n);gk=force_indx(3,n)
      if (      gi>=ngx1 .and. gi<=ngx2                          &
          .and. gj>=ngy1 .and. gj<=ngy2                          &
          .and. gk>=ngz1 .and. gk<=ngz2 ) nfrc=nfrc+1
   end do
   nmom=0
   do n=1,num_moment
      gi=moment_indx(1,n);gj=moment_indx(2,n);gk=moment_indx(3,n)
      if (      gi>=ngx1 .and. gi<=ngx2                          &
          .and. gj>=ngy1 .and. gj<=ngy2                          &
          .and. gk>=ngz1 .and. gk<=ngz2 ) nmom=nmom+1
   end do

   filenm=get_fnm_srcnode(n_i,n_j,n_k)
   call srcnode_skel(filenm,nfrc,ntwin_force,nmom,ntwin_moment)
   call nfseis_open(filenm,ncid)
  
   !force
   if (nfrc>0) then
      call nfseis_inq_varid(ncid,'force_indx',indxid)
      call nfseis_inq_varid(ncid,'force_axis',axisid)
      call nfseis_inq_varid(ncid,'force_start_time',t0id)
      call nfseis_inq_varid(ncid,'force_stf_time',stftid)
      call nfseis_inq_varid(ncid,'force_stf_freq',stffid)
      call nfseis_inq_varid(ncid,'Fx',fxid)
      call nfseis_inq_varid(ncid,'Fy',fyid)
      call nfseis_inq_varid(ncid,'Fz',fzid)
      call nfseis_put(ncid,stftid,frcstf_time,(/1/),(/ntwin_force/),(/1/))
      call nfseis_put(ncid,stffid,frcstf_freq,(/1/),(/ntwin_force/),(/1/))
      m=0
      do n=1,num_force
         gi=force_indx(1,n);gj=force_indx(2,n);gk=force_indx(3,n)
         if (      gi>=ngx1 .and. gi<=ngx2                          &
             .and. gj>=ngy1 .and. gj<=ngy2                          &
             .and. gk>=ngz1 .and. gk<=ngz2 ) then
            m=m+1
            call nfseis_put(ncid,fxid,ForceX(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
            call nfseis_put(ncid,fyid,ForceY(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
            call nfseis_put(ncid,fzid,ForceZ(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
            call nfseis_put(ncid,axisid,force_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
            call nfseis_put(ncid,t0id,force_t0(n),(/m/),(/1/),(/1/))
            call nfseis_put(ncid,indxid,           &
                 (/ out_i(swmpi_locli(gi,n_i)),    &
                    out_j(swmpi_loclj(gj,n_j)),    &
                    out_k(swmpi_loclk(gk,n_k)) /), &
                 (/1,m/),(/SEIS_GEO,1/),(/1,1/))
         end if
      end do
   end if
   !moment
   if (nmom>0) then
      call nfseis_inq_varid(ncid,'moment_indx',indxid)
      call nfseis_inq_varid(ncid,'moment_axis',axisid)
      call nfseis_inq_varid(ncid,'moment_start_time',t0id)
      call nfseis_inq_varid(ncid,'moment_stf_time',stftid)
      call nfseis_inq_varid(ncid,'moment_stf_freq',stffid)
      call nfseis_inq_varid(ncid,'Mxx',mxxid)
      call nfseis_inq_varid(ncid,'Myy',myyid)
      call nfseis_inq_varid(ncid,'Mzz',mzzid)
      call nfseis_inq_varid(ncid,'Mxy',mxyid)
      call nfseis_inq_varid(ncid,'Mxz',mxzid)
      call nfseis_inq_varid(ncid,'Myz',myzid)
      call nfseis_put(ncid,stftid,momstf_time,(/1/),(/ntwin_moment/),(/1/))
      call nfseis_put(ncid,stffid,momstf_freq,(/1/),(/ntwin_moment/),(/1/))
      m=0
      do n=1,num_moment
         gi=moment_indx(1,n);gj=moment_indx(2,n);gk=moment_indx(3,n)
         if (      gi>=ngx1 .and. gi<=ngx2                          &
             .and. gj>=ngy1 .and. gj<=ngy2                          &
             .and. gk>=ngz1 .and. gk<=ngz2 ) then
            m=m+1
            call nfseis_put(ncid,mxxid,MomTxx(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,myyid,MomTyy(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,mzzid,MomTzz(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,mxyid,MomTxy(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,mxzid,MomTxz(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,myzid,MomTyz(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
            call nfseis_put(ncid,axisid,moment_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
            call nfseis_put(ncid,t0id,moment_t0(n),(/m/),(/1/),(/1/))
            call nfseis_put(ncid,indxid,           &
                 (/ out_i(swmpi_locli(gi,n_i)),    &
                    out_j(swmpi_loclj(gj,n_j)),    &
                    out_k(swmpi_loclk(gk,n_k)) /), &
                 (/1,m/),(/SEIS_GEO,1/),(/1,1/))
         end if
      end do
   end if
   call nfseis_close(ncid)
end do
end do
end do

deallocate(gx)
deallocate(gy)
deallocate(gz)
call src_destroy

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

subroutine read_src_para(fnm_conf)
character (len=*),intent(in) :: fnm_conf
character (len=SEIS_STRLEN) :: mommech,stf_type,str
real(DP) :: strike,dip,rake
real(SP) :: d2m,f0,m0
integer fid,n,m

fid=1001
open(fid,file=trim(fnm_conf),status="old")

call string_conf(fid,1,'distance2meter',2,d2m)
call string_conf(fid,1,'src_hyper_height',2,src_hyper_height)
src_hyper_height=src_hyper_height*d2m

!force
call string_conf(fid,1,"number_of_force_source",2,num_force)
if (num_force>=1) then
   call string_conf(fid,1,"force_stf_window",2,ntwin_force)
   call string_conf(fid,1,"force_stf_type",2,stf_type)
   frcstf_id=stf_name2id(trim(stf_type))
   call src_alloc_force(num_force,ntwin_force)
   do m=1,ntwin_force
      call string_conf(fid,1,"force_stf_timefactor",m+1,frcstf_time(m))
      call string_conf(fid,1,"force_stf_freqfactor",m+1,frcstf_freq(m))
   end do
   call string_conf(fid,1,"<Force",2,str)
   do n=1,num_force
   do m=1,ntwin_force
      read(fid,*) force_axis(:,n),force_t0(n), &
          f0,ForceX(m,n),ForceY(m,n),ForceZ(m,n)
      ForceX(m,n)=ForceX(m,n)*f0
      ForceY(m,n)=ForceY(m,n)*f0
      ForceZ(m,n)=ForceZ(m,n)*f0
   end do
      force_axis(1:2,n)=force_axis(1:2,n)*PI/180.0_SP
      force_axis(3,n)=force_axis(3,n)*d2m
   end do
end if
!moment
call string_conf(fid,1,"number_of_moment_source",2,num_moment)
if (num_moment>=1) then
   call string_conf(fid,1,"moment_stf_window",2,ntwin_moment)
   call string_conf(fid,1,"moment_stf_type",2,stf_type)
   momstf_id=stf_name2id(trim(stf_type))
   call src_alloc_moment(num_moment,ntwin_moment)
   do m=1,ntwin_moment
      call string_conf(fid,1,"moment_stf_timefactor",m+1,momstf_time(m))
      call string_conf(fid,1,"moment_stf_freqfactor",m+1,momstf_freq(m))
   end do
   call string_conf(fid,1,"moment_mech_input",2,mommech)
   call string_conf(fid,1,"<Moment",2,str)
   if (trim(mommech)=='moment') then
      do n=1,num_moment
      do m=1,ntwin_moment
         read(fid,*) moment_axis(:,n),moment_t0(n),m0, &
           MomTxx(m,n),MomTyy(m,n),MomTzz(m,n),MomTxy(m,n),MomTxz(m,n),MomTyz(m,n)
         MomTxx(m,n)=m0*MomTxx(m,n);MomTyy(m,n)=m0*MomTyy(m,n);MomTzz(m,n)=m0*MomTzz(m,n)
         MomTxy(m,n)=m0*MomTxy(m,n);MomTxz(m,n)=m0*MomTxz(m,n);MomTyz(m,n)=m0*MomTyz(m,n)
      end do
         moment_axis(1:2,n)=moment_axis(1:2,n)*PI/180.0_SP
         moment_axis(3,n)=moment_axis(3,n)*d2m
      end do
   else
      do n=1,num_moment
      do m=1,ntwin_moment
         read(fid,*) moment_axis(:,n),moment_t0(n),m0,strike,dip,rake
         call angle2moment(strike,dip,rake, &
              MomTxx(m,n),MomTyy(m,n),MomTzz(m,n),MomTxy(m,n),MomTxz(m,n),MomTyz(m,n))
         MomTxx(m,n)=m0*MomTxx(m,n);MomTyy(m,n)=m0*MomTyy(m,n);MomTzz(m,n)=m0*MomTzz(m,n)
         MomTxy(m,n)=m0*MomTxy(m,n);MomTxz(m,n)=m0*MomTxz(m,n);MomTyz(m,n)=m0*MomTyz(m,n)
      end do
         moment_axis(1:2,n)=moment_axis(1:2,n)*PI/180.0_SP
         moment_axis(3,n)=moment_axis(3,n)*d2m
      end do
   end if
end if
close(fid)
end subroutine read_src_para

subroutine angle2moment(strike,dip,rake,Mxx,Myy,Mzz,Mxy,Mxz,Myz)
real(DP),intent(in) :: strike,dip,rake
real(SP),intent(out) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

real(DP) :: strike_pi,dip_pi,rake_pi
real(SP) :: M11,M22,M33,M12,M13,M23

dip_pi=dip/180.0_DP*PI
strike_pi=strike/180.0_DP*PI
rake_pi=rake/180.0_DP*PI
!in Aki and Richard's
M11=-(sin(dip_pi)*cos(rake_pi)*sin(2.0_DP*strike_pi)   &
         +sin(2.0_DP*dip_pi)*sin(rake_pi)*sin(strike_pi)**2)
M22=sin(dip_pi)*cos(rake_pi)*sin(2.0_DP*strike_pi)     &
         -sin(2.0_DP*dip_pi)*sin(rake_pi)*cos(strike_pi)**2
!Mzz=sin(2.0*dip_pi)*sin(rake_pi)
M33=-(Mxx+Myy)
M12=sin(dip_pi)*cos(rake_pi)*cos(2.0_DP*strike_pi)     &
         +0.5_DP*sin(2.0_DP*dip_pi)*sin(rake_pi)*sin(2.0_DP*strike_pi)
M13=-(cos(dip_pi)*cos(rake_pi)*cos(strike_pi)       &
         +cos(2.0_DP*dip_pi)*sin(rake_pi)*sin(strike_pi))
M23=-(cos(dip_pi)*cos(rake_pi)*sin(strike_pi)       &
         -cos(2.0_DP*dip_pi)*sin(rake_pi)*cos(strike_pi))
!to spherical
Mxx= M11; Myy=M22; Mzz=M33
Mxy=-M12; Mxz=M13; Myz=-M23
end subroutine angle2moment

!----------------------------------------------------------------

subroutine srcnode_skel(filenm,nfrc,ntwfrc,nmom,ntwmom)
character (len=*),intent(in) :: filenm
integer,intent(in) :: nfrc,ntwfrc,nmom,ntwmom

integer :: ncid,ierr,oldMode
integer :: geoid,vid,nmomid,ntwmomid,nfrcid,ntwfrcid

ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
     call nfseis_except(ierr,'srcnode_init:'//trim(filenm))
ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
     call nfseis_except(ierr,'set_fill in srcnode_skel')
! -- define dim
ierr=nf90_def_dim(ncid,'geo_dimension',SEIS_GEO,geoid)
     call nfseis_except(ierr,'geodim dim in srcnode_skel')
ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_force",nfrc)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_moment",nmom)
!force
if (nfrc>0) then
ierr=nf90_def_dim(ncid,'force_point',nfrc,nfrcid)
     call nfseis_except(ierr,'force_point dim in srcnode_skel')
ierr=nf90_def_dim(ncid,'force_time_window' ,ntwfrc,ntwfrcid)
     call nfseis_except(ierr,'force_time_window dim in srcnode_skel')
! -- define variable
ierr=nf90_def_var(ncid,'Fx', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'Fy', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'Fz', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_indx', nf90_int,(/ geoid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_axis', SEISNC_DATATYPE,(/ geoid, nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_start_time', SEISNC_DATATYPE,(/ nfrcid /),vid )
ierr=nf90_def_var(ncid,'force_stf_time', SEISNC_DATATYPE,(/ ntwfrcid /),vid )
ierr=nf90_def_var(ncid,'force_stf_freq', SEISNC_DATATYPE,(/ ntwfrcid /),vid )
ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_id",frcstf_id)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_type", &
      trim(stf_id2name(frcstf_id)) )
end if
!moment
if (nmom>0) then
ierr=nf90_def_dim(ncid,'moment_point',nmom,nmomid)
     call nfseis_except(ierr,'moment_point dim in srcnode_skel')
ierr=nf90_def_dim(ncid,'moment_time_window' ,ntwmom,ntwmomid)
     call nfseis_except(ierr,'moment_time_window dim in srcnode_skel')
! -- define variable
ierr=nf90_def_var(ncid,'Mxx',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Myy',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Mzz',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Mxy',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Mxz',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'Myz',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
ierr=nf90_def_var(ncid,'moment_indx', nf90_int,(/ geoid, nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_axis', SEISNC_DATATYPE,(/ geoid, nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_start_time', SEISNC_DATATYPE,(/ nmomid /),vid )
ierr=nf90_def_var(ncid,'moment_stf_time', SEISNC_DATATYPE,(/ ntwmomid /),vid )
ierr=nf90_def_var(ncid,'moment_stf_freq', SEISNC_DATATYPE,(/ ntwmomid /),vid )
ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_stf_id",momstf_id)
ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_stf_type", &
      trim(stf_id2name(momstf_id)) )
end if
ierr=nf90_enddef(ncid)
     call nfseis_except(ierr,'enddef in srcnode_skel')
ierr=nf90_close(ncid)
     call nfseis_except(ierr,'file close in srcnode_skel')
end subroutine srcnode_skel

end program seis3d_source

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
