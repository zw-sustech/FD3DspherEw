!******************************************************************************!
!*  This program calculate GSDF kernels for tau_p and tau_q.                  *!
!*  Following the code of Li Z, Po C, Zhang ZG and Shang Y                    *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

#define TomoObsConvSrc
!#define TomoFx
!#define TomoFy
#define TomoFz

!-----------------------------------------------------------------------------
program tomo_kernel
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
use io_mod

implicit none

INTEGER, PARAMETER :: TP = KIND(1.0D0)

character (len=SEIS_STRLEN) ::                 &
    fnm_main_conf,                             &
    fnm_tomo_conf,                             &
    pnm_obsinfo,pnm_obs,                       &
    pnm_result_source,                         &
    pnm_result_fx,pnm_result_fy,pnm_result_fz, &
    pnm_ker,fnm_prefix,                        &
    fnm_ker_log
integer ker_nt,ker_num,id_in_obs,nt1,nt2
integer id_of_obs, indx_in_obs, obs_nts,obs_ntt
real(SP) :: ker_dt,ker_win(2)
integer,allocatable :: ker_id(:)

real(SP),dimension(:),allocatable :: &
    Vx,Vy,Vz,Ux,Uy,Uz,T,S,           &
    Txx,Tyy,Tzz,Txy,Txz,Tyz,         &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,   &
    ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,   &
    Er,Es,Econ
real(kind=TP),dimension(:),allocatable :: Ka,Kb
real(SP) :: Kap,Kaq,Kbp,Kbq
real(SP) :: Vx0,Vy0,Vz0,Ux0,Uy0,Uz0
real(SP) :: src_m0
integer :: ker_nsift

character (len=SEIS_STRLEN) :: filenm
integer n_i,n_j,n_k
integer si,sj,sk,i,j,k,m,nid
logical iflag
integer,dimension(SEIS_GEO) :: &
   subs,subc,subt,sube,        &
   gsubs,gsubc,gsubt,gsube
integer p(1)

integer ncid
#ifdef TomoFx
integer kapxid,kaqxid,kbpxid,kbqxid
#endif
#ifdef TomoFy
integer kapyid,kaqyid,kbpyid,kbqyid
#endif
#ifdef TomoFz
integer kapzid,kaqzid,kbpzid,kbqzid
#endif

!----------------------------------------------------------------------!

call get_conf_name(fnm_conf)

! read kernel conf
fnm_tomo_conf='TomoKernel.conf'
call read_tomo_conf(fnm_tomo_conf)
fnm_conf=fnm_main_conf

call alloc_local(ker_nt)

call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)
call media_fnm_init(fnm_conf)
!call src_fnm_init(fnm_conf)

! alloc module
call grid_alloc
call media_alloc

! io
call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

nid=ker_id(1); ker_dt=snap_tinv(nid)*stept
T=(/0:ker_nt/)*ker_dt

obs_nts=snap_tinv(nid)/pt_tinv
obs_ntt=obs_nts

p=minloc(T,T>=ker_win(1)); nt1=p(1)-1
p=maxloc(T,T<=ker_win(2)); nt2=p(1)-1

! seismo on receiver
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
if (seismo_on_this(pnm_obsinfo,id_of_obs,indx_in_obs, &
                   n_i,n_j,n_k,m)) then
   filenm=get_fnm_seismo(pnm_obs,n_i,n_j,n_k)
   call retrieve_recvline(filenm,m,Vx(1:ker_nt),Vy(1:ker_nt),Vz(1:ker_nt), &
        obs_nts,ker_nt,obs_ntt)
   Vx(0)=0.0; Vy(0)=0.0; Vz(0)=0.0
   exit
end if
end do
end do
end do

#ifdef TomoObsConvSrc
call gen_stf(S,ker_nt,ker_dt,ker_nsift)
call cal_convolv(ker_nt,S,Vx,Ux,ker_dt,0,ker_nt); Vx=Ux
call cal_convolv(ker_nt,S,Vy,Uy,ker_dt,0,ker_nt); Vy=Uy
call cal_convolv(ker_nt,S,Vz,Uz,ker_dt,0,ker_nt); Vz=Uz
nt1=nt1+ker_nsift
nt2=nt2+ker_nsift
if (nt2>nt*stept/ker_dt) then
   print *, 'nt2 too large after src conv'
   print *, 'ndim, nt2=',nt*stept/ker_dt,nt2
   print *, 'src shift=',ker_nsift
   stop 1
end if
#endif

call vel2disp(Vx,Ux,ker_nt,ker_dt)
call vel2disp(Vy,Uy,ker_nt,ker_dt)
call vel2disp(Vz,Uz,ker_nt,ker_dt)

! normalization
Vx0=(0.5*Vx(nt1)**2+0.5*Vx(nt2)**2)*ker_dt
Vy0=(0.5*Vy(nt1)**2+0.5*Vy(nt2)**2)*ker_dt
Vz0=(0.5*Vz(nt1)**2+0.5*Vz(nt2)**2)*ker_dt
Ux0=(0.5*Ux(nt1)**2+0.5*Ux(nt2)**2)*ker_dt
Uy0=(0.5*Uy(nt1)**2+0.5*Uy(nt2)**2)*ker_dt
Uz0=(0.5*Uz(nt1)**2+0.5*Uz(nt2)**2)*ker_dt
do m=nt1+1,nt2-1
   Vx0=Vx0+Vx(m)**2.0*ker_dt
   Vy0=Vy0+Vy(m)**2.0*ker_dt
   Vz0=Vz0+Vz(m)**2.0*ker_dt
   Ux0=Ux0+Ux(m)**2.0*ker_dt
   Uy0=Uy0+Uy(m)**2.0*ker_dt
   Uz0=Uz0+Uz(m)**2.0*ker_dt
end do
Vx=Vx/Vx0/src_m0;Vy=Vy/Vy0/src_m0;Vz=Vz/Vz0/src_m0
Ux=Ux/Ux0/src_m0;Uy=Uy/Uy0/src_m0;Uz=Uz/Uz0/src_m0

! cal kernel
do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call media_import
do m=1,ker_num
   nid=ker_id(m)

   write(*,*) '       snap_id=',nid

   ! locate snap output in pnm_result_source
   filenm=get_fnm_snapnode_n( &
          trim(pnm_result_source),'norm_',nid,1,n_i,n_j,n_k)
   inquire(file=trim(filenm),exist=iflag)
   if (.not. iflag) cycle

   call nfseis_attget(filenm,'subs',subs)
   call nfseis_attget(filenm,'subc',subc)
   call nfseis_attget(filenm,'subt',subt)
   call nfseis_attget(filenm,'sube',sube)
   call nfseis_attget(filenm,'gsubs',gsubs)
   call nfseis_attget(filenm,'gsubc',gsubc)
   call nfseis_attget(filenm,'gsubt',gsubt)
   call nfseis_attget(filenm,'gsube',gsube)

   ! create kernel nc file
   filenm=get_fnm_snapnode_n(pnm_ker,'kernel_',nid,0,n_i,n_j,n_k)
   call nfseis_data_create(filenm,subc(1),subc(2),subc(3), &
        "Finite frequency kernel")
   call nfseis_attput(filenm,'subs',subs)
   call nfseis_attput(filenm,'subc',subc)
   call nfseis_attput(filenm,'subt',subt)
   call nfseis_attput(filenm,'sube',sube)
   call nfseis_attput(filenm,'gsubs',gsubs)
   call nfseis_attput(filenm,'gsubc',gsubc)
   call nfseis_attput(filenm,'gsubt',gsubt)
   call nfseis_attput(filenm,'gsube',gsube)
#ifdef TomoFx
   call nfseis_data_addvar(filenm,'Kapx')
   call nfseis_data_addvar(filenm,'Kaqx')
   call nfseis_data_addvar(filenm,'Kbpx')
   call nfseis_data_addvar(filenm,'Kbqx')
#endif
#ifdef TomoFy
   call nfseis_data_addvar(filenm,'Kapy')
   call nfseis_data_addvar(filenm,'Kaqy')
   call nfseis_data_addvar(filenm,'Kbpy')
   call nfseis_data_addvar(filenm,'Kbqy')
#endif
#ifdef TomoFz
   call nfseis_data_addvar(filenm,'Kapz')
   call nfseis_data_addvar(filenm,'Kaqz')
   call nfseis_data_addvar(filenm,'Kbpz')
   call nfseis_data_addvar(filenm,'Kbqz')
#endif

   call nfseis_open(filenm,ncid)

#ifdef TomoFx
   call nfseis_inq_varid(ncid,'Kapx',kapxid)
   call nfseis_inq_varid(ncid,'Kaqx',kaqxid)
   call nfseis_inq_varid(ncid,'Kbpx',kbpxid)
   call nfseis_inq_varid(ncid,'Kbqx',kbqxid)
#endif
#ifdef TomoFy
   call nfseis_inq_varid(ncid,'Kapy',kapyid)
   call nfseis_inq_varid(ncid,'Kaqy',kaqyid)
   call nfseis_inq_varid(ncid,'Kbpy',kbpyid)
   call nfseis_inq_varid(ncid,'Kbqy',kbqyid)
#endif
#ifdef TomoFz
   call nfseis_inq_varid(ncid,'Kapz',kapzid)
   call nfseis_inq_varid(ncid,'Kaqz',kaqzid)
   call nfseis_inq_varid(ncid,'Kbpz',kbpzid)
   call nfseis_inq_varid(ncid,'Kbqz',kbqzid)
#endif

   subs(1)=inn_i(subs(1));subs(2)=inn_j(subs(2));subs(3)=inn_k(subs(3))
   sube(1)=inn_i(sube(1));sube(2)=inn_j(sube(2));sube(3)=inn_k(sube(3))

sk=0
do k=subs(3),sube(3),subt(3)
   sk=sk+1; sj=0
   write(*,*) '          k=',k
do j=subs(2),sube(2),subt(2)
   sj=sj+1; si=0
do i=subs(1),sube(1),subt(1)
   si=si+1

   ! read inner source stress
   fnm_prefix=trim(pnm_result_source)                            &
       //'/'//trim(io_enum('norm_',nid))                         &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txx',Txx(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyy',Tyy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tzz',Tzz(1:nt2),nt2)
   Txx(0)=0.0;Tyy(0)=0.0;Tzz(0)=0.0;
   fnm_prefix=trim(pnm_result_source)                            &
       //'/'//trim(io_enum('shear_',nid))                        &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txy',Txy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txz',Txz(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyz',Tyz(1:nt2),nt2)
   Txy(0)=0.0;Txz(0)=0.0;Tyz(0)=0.0;
   call stress2strain(Txx,Tyy,Tzz,Txy,Txz,Tyz,                   &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,                           &
        ker_nt, nt2, lambda(i,j,k),mu(i,j,k))

#ifdef TomoFx
   ! read green stress caused by Fx
   fnm_prefix=trim(pnm_result_fx)                                &
       //'/'//trim(io_enum('norm_',nid))                         &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txx',Txx(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyy',Tyy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tzz',Tzz(1:nt2),nt2)
   Txx(0)=0.0;Tyy(0)=0.0;Tzz(0)=0.0;
   fnm_prefix=trim(pnm_result_fx)                                &
       //'/'//trim(io_enum('shear_',nid))                        &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txy',Txy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txz',Txz(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyz',Tyz(1:nt2),nt2)
   Txy(0)=0.0;Txz(0)=0.0;Tyz(0)=0.0;
#ifndef TomoObsConvSrc
   call corr_green(Txx,Tyy,Tzz,Txy,Txz,Tyz,stf_t0,ker_dt,ker_nt)
#endif
   call stress2strain(Txx,Tyy,Tzz,Txy,Txz,Tyz,     &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,             &
        ker_nt, nt2,lambda(i,j,k),mu(i,j,k))
   ! cal
   call cal_kernel(Vx,Ux,                          &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,             &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,             &
        lambda(i,j,k),mu(i,j,k),                   &
        ker_nt,nt1,nt2,ker_dt,Kap,Kaq,Kbp,Kbq)
   !put
   call nfseis_put(ncid,kapxid,Kap,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kaqxid,Kaq,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kbpxid,Kbp,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kbqxid,Kbq,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
#endif

#ifdef TomoFy
   ! read green stress caused by Fy
   fnm_prefix=trim(pnm_result_fy)                                &
       //'/'//trim(io_enum('norm_',nid))                         &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txx',Txx(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyy',Tyy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tzz',Tzz(1:nt2),nt2)
   Txx(0)=0.0;Tyy(0)=0.0;Tzz(0)=0.0;
   fnm_prefix=trim(pnm_result_fy)                                &
       //'/'//trim(io_enum('shear_',nid))                        &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txy',Txy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txz',Txz(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyz',Tyz(1:nt2),nt2)
   Txy(0)=0.0;Txz(0)=0.0;Tyz(0)=0.0;
#ifndef TomoObsConvSrc
   call corr_green(Txx,Tyy,Tzz,Txy,Txz,Tyz,stf_t0,ker_dt,ker_nt)
#endif
   call stress2strain(Txx,Tyy,Tzz,Txy,Txz,Tyz,     &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,             &
        ker_nt, nt2,lambda(i,j,k),mu(i,j,k))
   ! cal
   call cal_kernel(Vy,Uy,                          &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,             &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,             &
        lambda(i,j,k),mu(i,j,k),                   &
        ker_nt,nt1,nt2,ker_dt,Kap,Kaq,Kbp,Kbq)
   !put
   call nfseis_put(ncid,kapyid,Kap,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kaqyid,Kaq,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kbpyid,Kbp,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kbqyid,Kbq,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
#endif

#ifdef TomoFz
   ! read green stress caused by Fz
   fnm_prefix=trim(pnm_result_fz)                                &
       //'/'//trim(io_enum('norm_',nid))                         &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txx',Txx(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyy',Tyy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tzz',Tzz(1:nt2),nt2)
   Txx(0)=0.0;Tyy(0)=0.0;Tzz(0)=0.0;
   fnm_prefix=trim(pnm_result_fz)                                &
       //'/'//trim(io_enum('shear_',nid))                        &
       //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txy',Txy(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Txz',Txz(1:nt2),nt2)
   call retrieve_snap_seis(fnm_prefix,si,sj,sk,'Tyz',Tyz(1:nt2),nt2)
   Txy(0)=0.0;Txz(0)=0.0;Tyz(0)=0.0;
#ifndef TomoObsConvSrc
   call corr_green(Txx,Tyy,Tzz,Txy,Txz,Tyz,stf_t0,ker_dt,ker_nt)
#endif
   call stress2strain(Txx,Tyy,Tzz,Txy,Txz,Tyz,     &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,             &
        ker_nt, nt2,lambda(i,j,k),mu(i,j,k))
   ! cal
   call cal_kernel(Vz,Uz,                          &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,             &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,             &
        lambda(i,j,k),mu(i,j,k),                   &
        ker_nt,nt1,nt2,ker_dt,Kap,Kaq,Kbp,Kbq)
   !put
   call nfseis_put(ncid,kapzid,Kap,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kaqzid,Kaq,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kbpzid,Kbp,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
   call nfseis_put(ncid,kbqzid,Kbq,(/si,sj,sk/),(/1,1,1/),(/1,1,1/))
#endif

end do !i
end do !j
end do !k

   call nfseis_close(ncid)

end do !m

end do !n_k
end do
end do !n_i

call dealloc_local

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

subroutine alloc_local(nt)
  integer nt
  allocate(Vx(0:nt)); Vx=0.0
  allocate(Vy(0:nt)); Vy=0.0
  allocate(Vz(0:nt)); Vz=0.0
  allocate(Ux(0:nt)); Ux=0.0
  allocate(Uy(0:nt)); Uy=0.0
  allocate(Uz(0:nt)); Uz=0.0
  allocate(T (0:nt)); T =0.0
  allocate(S (0:nt)); S =0.0
  allocate(Txx(0:nt),Tyy(0:nt),Tzz(0:nt),Txy(0:nt),Txz(0:nt),Tyz(0:nt))
  allocate(ExxS(0:nt),EyyS(0:nt),EzzS(0:nt),ExyS(0:nt),ExzS(0:nt),EyzS(0:nt))
  allocate(ExxR(0:nt),EyyR(0:nt),EzzR(0:nt),ExyR(0:nt),ExzR(0:nt),EyzR(0:nt))
  allocate(Es(0:nt),Er(0:nt),Econ(0:nt))
  allocate(Ka(0:nt),Kb(0:nt))
end subroutine alloc_local

subroutine dealloc_local
  if (allocated(Vx)) deallocate(Vx)
  if (allocated(Vy)) deallocate(Vy)
  if (allocated(Vz)) deallocate(Vz)
  if (allocated(Ux)) deallocate(Ux)
  if (allocated(Uy)) deallocate(Uy)
  if (allocated(Uz)) deallocate(Uz)
  if (allocated(T )) deallocate(T )
  if (allocated(Txx)) deallocate(Txx)
  if (allocated(Tyy)) deallocate(Tyy)
  if (allocated(Tzz)) deallocate(Tzz)
  if (allocated(Txy)) deallocate(Txy)
  if (allocated(Txz)) deallocate(Txz)
  if (allocated(Tyz)) deallocate(Tyz)
  if (allocated(ExxS)) deallocate(ExxS)
  if (allocated(EyyS)) deallocate(EyyS)
  if (allocated(EzzS)) deallocate(EzzS)
  if (allocated(ExyS)) deallocate(ExyS)
  if (allocated(ExzS)) deallocate(ExzS)
  if (allocated(EyzS)) deallocate(EyzS)
  if (allocated(ExxR)) deallocate(ExxR)
  if (allocated(EyyR)) deallocate(EyyR)
  if (allocated(EzzR)) deallocate(EzzR)
  if (allocated(ExyR)) deallocate(ExyR)
  if (allocated(ExzR)) deallocate(ExzR)
  if (allocated(EyzR)) deallocate(EyzR)
  if (allocated(Ka)) deallocate(Ka)
  if (allocated(Kb)) deallocate(Kb)
end subroutine dealloc_local

subroutine read_tomo_conf(fnm_conf)
character (len=*) :: fnm_conf
character (len=SEIS_STRLEN) :: fnm_src
character (len=SEIS_STRLEN) :: stf_type
integer fid,n
fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'fnm_main_conf',2,fnm_main_conf)
call string_conf(fid,1,'fnm_ker_log',2,fnm_ker_log)

call string_conf(fid,1,'pnm_obsinfo',2,pnm_obsinfo)
call string_conf(fid,1,'pnm_obs',2,pnm_obs)
call string_conf(fid,1,'id_of_obs',2,id_of_obs)
call string_conf(fid,1,'indx_in_obs',2,indx_in_obs)
call string_conf(fid,1,'time_window',2,ker_win(1))
call string_conf(fid,1,'time_window',3,ker_win(2))
call string_conf(fid,1,'fnm_src_conf',2,fnm_src)
call string_conf(fid,1,"m0",2,src_m0)

call string_conf(fid,1,'pnm_result_source',2,pnm_result_source)
call string_conf(fid,1,'pnm_result_fx',2,pnm_result_fx)
call string_conf(fid,1,'pnm_result_fy',2,pnm_result_fy)
call string_conf(fid,1,'pnm_result_fz',2,pnm_result_fz)
call string_conf(fid,1,'pnm_tomo_kernel',2,pnm_ker)
call string_conf(fid,1,'kernel_snap_nt',2,ker_nt)
call string_conf(fid,1,'num_kernel',2,ker_num)
allocate(ker_id(ker_num))
do n=1,ker_num
   call string_conf(fid,1,'kernel_snap_id',n+1,ker_id(n))
end do
close(fid)

open(fid,file=trim(fnm_src),status='old')
call string_conf(fid,1,"ntwin",2,ntwin)
allocate(twin(2,ntwin))
do n=1,ntwin
   call string_conf(fid,1,"twin",2*n,twin(1,n))
   call string_conf(fid,1,"twin",2*n+1,twin(2,n))
end do
call string_conf(fid,1,"stf",2,stf_type)
select case (trim(stf_type))
case ('bell_int')
     stf_t0=twin(2,1)-twin(1,1)
     stf_alpha=0.0
     flag_stf_type=SIG_SVF_BELL
case ('bell')
     stf_t0=twin(2,1)-twin(1,1)
     stf_alpha=0.0
     flag_stf_type=SIG_STF_BELL
case ('triangle_int')
     stf_t0=twin(2,1)-twin(1,1)
     stf_alpha=0.0
     flag_stf_type=SIG_SVF_TRIANGLE
case ('ricker')
     call string_conf(fid,1,"stf_para",2,stf_alpha)
     call string_conf(fid,1,"stf_para",3,stf_t0)
     flag_stf_type=SIG_STF_RICKER
case ('B(shift)')
     call string_conf(fid,1,"stf_para",2,stf_t0)
     flag_stf_type=SIG_STF_BSHIFT
     print *, 'B(shift) not finished'
     stop 1
case default
     print *, "Have you told me how to generate the STF of ", trim(stf_type)
     stop 1
end select
close(fid)
end subroutine read_tomo_conf

subroutine cal_kernel_conv(V,U,           &
           ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
           ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
           lam,miu,                       &
           nt,nt1,nt2,dt,Kap,Kaq,Kbp,Kbq)

integer :: nt,nt1,nt2
real(SP),dimension(0:nt) :: V,U,    &
     ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
     ExxR,EyyR,EzzR,ExyR,ExzR,EyzR
real(SP) :: miu,lam
real(SP) :: dt,Kap,Kaq,Kbp,Kbq

integer n
real(SP) Vprho,Vsrho

!p=minloc(T,T>=t1); nt1=p(1)
!p=maxloc(T,T<=t2); nt2=p(1)

do n=0,nt2
   Er(n)=ExxR(n)+EyyR(n)+EzzR(n)
   Es(n)=ExxS(n)+EyyS(n)+EzzS(n)
end do

! for Kap
call cal_convolv(nt,Er,Es,Econ,dt,nt1,nt2)

! for Kbp
do n=nt1,nt2
   Es(n)= -Econ(n)
end do
call cal_convolv(nt,ExxR,ExxS,Er,dt,nt1,nt2)
do n=nt1,nt2
   Es(n)=Es(n)+Er(n)
end do
call cal_convolv(nt,EyyR,EyyS,Er,dt,nt1,nt2)
do n=nt1,nt2
   Es(n)=Es(n)+Er(n)
end do
call cal_convolv(nt,EzzR,EzzS,Er,dt,nt1,nt2)
do n=nt1,nt2
   Es(n)=Es(n)+Er(n)
end do
call cal_convolv(nt,ExyR,ExyS,Er,dt,nt1,nt2)
do n=nt1,nt2
   Es(n)=Es(n)+Er(n)*2.0
end do
call cal_convolv(nt,ExzR,ExzS,Er,dt,nt1,nt2)
do n=nt1,nt2
   Es(n)=Es(n)+Er(n)*2.0
end do
call cal_convolv(nt,EyzR,EyzS,Er,dt,nt1,nt2)
do n=nt1,nt2
   Es(n)=Es(n)+Er(n)*2.0
   Es(n)=2.0*Es(n)
end do
 
!   Forming GSDF kernels for tau_p and tau_q.
Kap=0.5*dt*(V(nt1)*Econ(nt1)+V(nt2)*Econ(nt2))
Kaq=0.5*dt*(U(nt1)*Econ(nt1)+U(nt2)*Econ(nt2))
Kbp=0.5*dt*(V(nt1)*Es(nt1)+V(nt2)*Es(nt2))
Kbq=0.5*dt*(U(nt1)*Es(nt1)+U(nt2)*Es(nt2))
do n=nt1+1,nt2-1
   Kap=Kap+V(n)*Econ(n)*dt
   Kaq=Kaq+U(n)*Econ(n)*dt
   Kbp=Kbp+V(n)*Es(n)*dt
   Kbq=Kbq+U(n)*Es(n)*dt
end do

!Vprho=sqrt( (lam+2.0*miu)*dens )
!Vsrho=sqrt( miu*dens )
Vprho= lam+2.0*miu
Vsrho= miu
Kap=Kap*2.0*Vprho
Kaq=-Kaq*2.0*Vprho
Kbp=Kbp*2.0*Vsrho
Kbq=-Kbq*2.0*Vsrho
end subroutine cal_kernel_conv

subroutine cal_kernel(V,U,                &
           ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
           ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
           lam,miu,                       &
           nt,nt1,nt2,dt,Kap,Kaq,Kbp,Kbq)

integer :: nt,nt1,nt2
real(SP),dimension(0:nt) :: V,U,          &
     ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,       &
     ExxR,EyyR,EzzR,ExyR,ExzR,EyzR
real(SP) :: miu,lam
real(SP) :: dt,Kap,Kaq,Kbp,Kbq

integer n,i
real(SP) Vprho,Vsrho
real(kind=TP) :: x1,x2

x1=0.0_TP
x2=0.0_TP

do n=nt1,nt2

x1= 0.5_TP*dt*(ExxR(0)+EyyR(0)+EzzR(0))*(ExxS(n)+EyyS(n)+EzzS(n)) &
   +0.5_TP*dt*(ExxR(n)+EyyR(n)+EzzR(n))*(ExxS(0)+EyyS(0)+EzzS(0))
x2= 0.5_TP*dt*( ExxR(0)*ExxS(n)+EyyR(0)*EyyS(n)+EzzR(0)*EzzS(n)                         &
               +2.0_TP*ExyR(0)*ExyS(n)+2.0_TP*ExzR(0)*ExzS(n)+2.0_TP*EyzR(0)*EyzS(n) )  &
   +0.5_TP*dt*( ExxR(n)*ExxS(0)+EyyR(n)*EyyS(0)+EzzR(n)*EzzS(0)                         &
               +2.0_TP*ExyR(n)*ExyS(0)+2.0_TP*ExzR(n)*ExzS(0)+2.0_TP*EyzR(n)*EyzS(0) )
do i=1,n-1
   x1=x1+1.0_TP*dt*(ExxR(i)+EyyR(i)+EzzR(i))*(ExxS(n-i)+EyyS(n-i)+EzzS(n-i))
   x2=x2+1.0_TP*dt*( ExxR(i)*ExxS(n-i)+EyyR(i)*EyyS(n-i)+EzzR(i)*EzzS(n-i)     &
               +2.0_TP*ExyR(i)*ExyS(n-i)+2.0_TP*ExzR(i)*ExzS(n-i)+2.0_TP*EyzR(i)*EyzS(n-i) )
end do
Ka(n)=x1
Kb(n)=2.0_TP*(x2-x1)
end do

!   Forming GSDF kernels for tau_p and tau_q.
Kap=0.5*dt*(V(nt1)*Ka(nt1)+V(nt2)*Ka(nt2))
Kaq=0.5*dt*(U(nt1)*Ka(nt1)+U(nt2)*Ka(nt2))
Kbp=0.5*dt*(V(nt1)*Kb(nt1)+V(nt2)*Kb(nt2))
Kbq=0.5*dt*(U(nt1)*Kb(nt1)+U(nt2)*Kb(nt2))
do n=nt1+1,nt2-1
   Kap=Kap+V(n)*Ka(n)*dt
   Kaq=Kaq+U(n)*Ka(n)*dt
   Kbp=Kbp+V(n)*Kb(n)*dt
   Kbq=Kbq+U(n)*Kb(n)*dt
end do

Vprho= lam+2.0*miu
Vsrho= miu
Kap=Kap*2.0*Vprho
Kaq=-Kaq*2.0*Vprho
Kbp=Kbp*2.0*Vsrho
Kbq=-Kbq*2.0*Vsrho
end subroutine cal_kernel

subroutine cal_convolv(nt,U,V,W,dt,nt1,nt2)
integer nt,nt1,nt2
real(SP),dimension(0:nt) :: U,V,W
real(SP) dt
integer i,j

!W=0.0
do i=nt1,nt2
   W(i)=0.5*dt*U(0)*V(i)+0.5*dt*U(i)*V(0)
   do j=1,i-1
      W(i)=W(i)+dt*U(j)*V(i-j)
   end do
end do
end subroutine cal_convolv

subroutine vel2disp(V,U,ntpt,dt)
integer ntpt
real(SP),dimension(0:ntpt) :: V,U
real(SP) dt
integer n

!U=(sum(V(1:nt-1))+V(nt)*0.5)*dt
U(0)=0.0
do n=1,ntpt
   U(n)=U(n-1)+0.5*(V(n-1)+V(n))*dt
end do
end subroutine vel2disp

subroutine corr_green(Txx,Tyy,Tzz,Txy,Txz,Tyz,t0,dt,num)
real(SP),dimension(:) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
real(SP) :: t0,dt
integer num
integer n0

n0=nint(t0/dt)
Txx(1:num-n0)=Txx(n0+1:num)
Tyy(1:num-n0)=Tyy(n0+1:num)
Tzz(1:num-n0)=Tzz(n0+1:num)
Txy(1:num-n0)=Txy(n0+1:num)
Txz(1:num-n0)=Txz(n0+1:num)
Tyz(1:num-n0)=Tyz(n0+1:num)
!Txx=Txx/src_m0; Tyy=Tyy/src_m0; Tzz=Tzz/src_m0
!Txy=Txy/src_m0; Txz=Txz/src_m0; Tyz=Tyz/src_m0
end subroutine corr_green

subroutine stress2strain(Txx,Tyy,Tzz,Txy,Txz,Tyz, &
   Exx,Eyy,Ezz,Exy,Exz,Eyz,                       &
   nt,nt2,lam,miu)
integer nt,nt2
real(SP),dimension(0:nt) :: Txx,Tyy,Tzz,Txy,Txz,Tyz
real(SP),dimension(0:nt) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
real(SP) lam,miu
real(SP) E1,E2,E3
integer n

E1=(lam+miu)/(miu*(3.0*lam+2.0*miu))
E2=-lam/(2.0*miu*(3.0*lam+2.0*miu))
E3=1.0/miu
Exx(0)=0.0;Eyy(0)=0.0;Ezz(0)=0.0
Exy(0)=0.0;Exz(0)=0.0;Eyz(0)=0.0
do n=1,nt2
   Exx(n)=E1*Txx(n)+E2*Tyy(n)+E2*Tzz(n)
   Eyy(n)=E2*Txx(n)+E1*Tyy(n)+E2*Tzz(n)
   Ezz(n)=E2*Txx(n)+E2*Tyy(n)+E1*Tzz(n)
   Exy(n)=0.5*E3*Txy(n)
   Exz(n)=0.5*E3*Txz(n)
   Eyz(n)=0.5*E3*Tyz(n)
end do
end subroutine stress2strain

subroutine gen_stf(S,nt,dt,n0)
integer nt,n0
real(SP),dimension(0:nt) :: S
real(SP) :: riset,rate,dt,t
integer n
S=0.0
riset=twin(2,1)-twin(1,1)
do n=1,nt
   t=n*dt

select case(flag_stf_type)
case (SIG_STF_BELL)
   rate= fun_bell(t,riset)
   n0=nint(riset/2.0/dt)
case (SIG_STF_RICKER)
   rate=fun_ricker(t,stf_alpha,stf_t0)
   n0=nint(stf_t0/dt)
case (SIG_STF_GAUSS)
   rate=fun_gauss(t,stf_alpha,stf_t0)
   n0=nint(stf_t0/dt)
case default
   print *, "Have you told me how to generate the STF of ", flag_stf_type
   stop 1
end select
if (abs(rate)<SEIS_ZERO) rate=0.0
  
   S(n)=rate

end do
print *, sum(S)*dt
end subroutine gen_stf

end program tomo_kernel
