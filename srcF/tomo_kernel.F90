!******************************************************************************!
!*  This program calculates finite-frequency kernels of tau_p and tau_q.      *!
!*  Following the code of L ZHAO, P CHEN, ZG ZHANG and Y SHEN                 *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
program tomo_kernel
!-----------------------------------------------------------------------------

#define VERBOSE
#define USEMPI

use constants_mod
use string_mod
use para_mod
use io_mod
use media_mod
use nfseis_mod
use src_mod
use mpi_mod
#ifdef USEMPI
use mpi
#endif

implicit none

character (len=SEIS_STRLEN) ::   &
    fnm_main_conf,fnm_tomo_conf, &
    pnm_wave,pnm_sgt,pnm_ker,    &
    fnm_record,fnm_ker,ker_prefix
character (len=SEIS_STRLEN) :: filenm

integer,dimension(SEIS_GEO) :: blklen,blknum
integer :: kid,knt,nt1,nt2
real(SP) :: kdt,kwin(2)

real(SP),dimension(:),allocatable :: V,U,T
real(SP) :: V0,U0,src_m0

real(SP),dimension(:,:,:,:),allocatable :: &
    ExxS,EyyS,EzzS,ExyS,ExzS,EyzS ,        &
    ExxR,EyyR,EzzR,ExyR,ExzR,EyzR

real(SP),dimension(:),allocatable :: Ka,Kb
real(SP),dimension(:),allocatable :: E11R,E22R,E33R,E12R,E13R,E23R
real(SP),dimension(:),allocatable :: E11S,E22S,E33S,E12S,E13S,E23S
real(SP),dimension(:,:,:),allocatable :: Kap,Kaq,Kbp,Kbq

logical :: flag_filter
integer :: NOD,NOD1,NEL
real(DP),dimension(:),allocatable :: filt_a,filt_b
real(SP),dimension(:),allocatable :: vecx,vecy,vecw

integer :: n_i,n_j,n_k
integer :: i,j,k,m,n,mt,ierr
integer,dimension(SEIS_GEO) :: &
   bsubs,bsubc,bsubt,          &
    subs, subc, subt
integer :: p(1)
#ifdef USEMPI
integer,dimension(MPI_STATUS_SIZE) :: istatus
#endif

integer :: ncid,ncidS,ncidF
integer :: kapid,kaqid,kbpid,kbqid
integer :: TxxSid,TyySid,TzzSid,TxySid,TxzSid,TyzSid
integer :: TxxFid,TyyFid,TzzFid,TxyFid,TxzFid,TyzFid

!----------------------------------------------------------------------

#ifdef USEMPI
call MPI_INIT(ierr)
#endif

call get_conf_name(fnm_conf)

! read kernel conf
fnm_tomo_conf='TomoKernel.conf'
call read_tomo_conf(fnm_tomo_conf)
fnm_conf=fnm_main_conf

call para_init(fnm_conf)
call swmpi_init(fnm_conf)

#ifdef MediaMPI
call swmpi_cart_creat
call swmpi_reinit_para
call swmpi_datatype
#else
call swmpi_set_gindx(0,0,0)
#endif

call media_fnm_init(fnm_conf)

call io_init(fnm_conf)
call io_snap_read(fnm_conf)
call io_pt_read(fnm_conf)

!----------------------------------------------------------------------
! seismo on receiver
!#ifdef USEMPI
!if (seismo_on_this(pnm_obsinfo,id_of_obs,indx_in_obs, &
!                   thisid(1),thisid(2),thisid(3),n)) then
!   filenm=get_fnm_seismo(pnm_obs,thisid(1),thisid(2),thisid(3))
!   i=snap_tinv(kid)/pt_tinv; ! recv type ouput may be denser than snap
!   call retrieve_recvline(filenm,n,varnm_obs,Vz,i,knt,i)
!   do n=0,dims(1)*dims(2)*dims(3)-1
!   if (n/=myid) then
!      call MPI_SEND(Vz,knt,MPI_REAL,n,3,MPI_COMM_WORLD,ierr)
!   end if
!   end do
!else
!   call MPI_RECV(Vz,knt,MPI_REAL,MPI_ANY_SOURCE,3,MPI_COMM_WORLD,istatus,ierr)
!end if
!#endif

kdt=snap_tinv(kid)*stept
if (abs(T(2)-T(1)-kdt)>SEIS_ZERO) call error_except('dt /= stept*tinv')

p=minloc(T,T>=kwin(1)); nt1=p(1)
p=maxloc(T,T<=kwin(2)); nt2=p(1)

! normalization
V0=(0.5*V(nt1)**2+0.5*V(nt2)**2)*kdt; U0=(0.5*U(nt1)**2+0.5*U(nt2)**2)*kdt
do m=nt1+1,nt2-1
   V0=V0+V(m)**2.0*kdt; U0=U0+U(m)**2.0*kdt
end do
V=V/V0/src_m0; U=U/U0/src_m0

!-----------------------------------------------------------------------------
#ifndef USEMPI
print *, 'input mpi id:'
read *, n_i,n_j,n_k
call swmpi_change_fnm(n_i,n_j,n_k)
call swmpi_set_gindx(n_i,n_j,n_k)
thisid=(/ n_i,n_j,n_k /)
#else
n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#endif

call io_snap_locate(n_i,n_j,n_k)

kernode_if : if (snap_ishere(kid) .and. sgt_out(kid)) then

call ker_para_reinit
call alloc_local(blklen(1),blklen(2),blklen(3),knt)
call alloc_media_local(blklen(1),blklen(2),blklen(3))

! create kernel nc file
fnm_ker=get_fnm_snapnode_n(pnm_ker,ker_prefix,kid,0,n_i,n_j,n_k)
call nfseis_grid3d_def(fnm_ker,                                &
     snap_subc(1,kid),snap_subc(2,kid), snap_subc(3,kid),ncid, &
     "Finite frequency kernel for record "//trim(fnm_record))
call nfseis_grid3d_defvar(ncid,'kernel_phase_Vp'    ,kapid)
call nfseis_grid3d_defvar(ncid,'kernel_amplitude_Vp',kaqid)
call nfseis_grid3d_defvar(ncid,'kernel_phase_Vs'    ,kbpid)
call nfseis_grid3d_defvar(ncid,'kernel_amplitude_Vs',kbqid)
call nfseis_grid3d_enddef(ncid)

! loop each block
do k=1,blknum(3)
do j=1,blknum(2)
do i=1,blknum(1)

#ifdef VERBOSE
   write(*,"(a7,3(i4,a1,i4.4),a9,3(i2.2))") &
       ' block:',                           &
       i,'/',blknum(1),                     &
       j,'/',blknum(2),                     &
       k,'/',blknum(3),                     &
       ', thisid=',n_i,n_j,n_k
#endif

   bsubs=(/ (i-1)*blklen(1)+1,(j-1)*blklen(2)+1,(k-1)*blklen(3)+1 /)
   bsubc=(/ min(blklen(1),snap_subc(1,kid)-bsubs(1)+1),  &
            min(blklen(2),snap_subc(2,kid)-bsubs(2)+1),  &
            min(blklen(3),snap_subc(3,kid)-bsubs(3)+1) /)
   bsubt=(/ 1,1,1 /)
   subc=bsubc; subt=snap_subt(:,kid)*bsubt
   subs=snap_subs(:,kid)+(bsubs-1)*subt

if (product(bsubc)*knt/=size(ExxS)) then
   call realloc_local(bsubc(1),bsubc(2),bsubc(3),knt)
end if

! load modeling field
if (snap_tcnt(kid)<knt) then
   m=0; n=0; mt=knt
   do 
      if (m>=knt) exit
      n=n+1
      filenm=get_fnm_snapnode_n(pnm_wave,'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidS)
      call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
      call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
      call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
      call nfseis_inq_varid(ncidS, 'Txy', TxySid)
      call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
      call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)

      filenm=get_fnm_snapnode_n(pnm_sgt,'sgt_',kid,n,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
      if (m+mt>knt) mt=knt-m

      ierr=nf90_get_var(ncidS,TxxSid,ExxS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyySid,EyyS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TzzSid,EzzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxySid,ExyS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TxzSid,ExzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidS,TyzSid,EyzS(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))

      ierr=nf90_get_var(ncidF,TxxFid,ExxR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyyFid,EyyR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TzzFid,EzzR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxyFid,ExyR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TxzFid,ExzR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      ierr=nf90_get_var(ncidF,TyzFid,EyzR(:,:,:,m+1:m+mt),(/bsubs,1/),(/bsubc,mt/),(/bsubt,1/))
      if (ierr /= nf90_noerr) &
      call nfseis_except(ierr,'tomo_kernel_mpi: get equavilent stress fail')
      m=m+mt
   end do
   call nfseis_close(ncidS)
   call nfseis_close(ncidF)
else
   if (i==1 .and. j==1 .and. k==1) then
      filenm=get_fnm_snapnode_n(pnm_wave,'sgt_',kid,1,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidS)
      call nfseis_inq_varid(ncidS, 'Txx', TxxSid)
      call nfseis_inq_varid(ncidS, 'Tyy', TyySid)
      call nfseis_inq_varid(ncidS, 'Tzz', TzzSid)
      call nfseis_inq_varid(ncidS, 'Txy', TxySid)
      call nfseis_inq_varid(ncidS, 'Txz', TxzSid)
      call nfseis_inq_varid(ncidS, 'Tyz', TyzSid)
   
      filenm=get_fnm_snapnode_n(pnm_sgt,'sgt_',kid,1,n_i,n_j,n_k)
      call nfseis_open(filenm,ncidF)
      call nfseis_inq_varid(ncidF, 'Txx', TxxFid)
      call nfseis_inq_varid(ncidF, 'Tyy', TyyFid)
      call nfseis_inq_varid(ncidF, 'Tzz', TzzFid)
      call nfseis_inq_varid(ncidF, 'Txy', TxyFid)
      call nfseis_inq_varid(ncidF, 'Txz', TxzFid)
      call nfseis_inq_varid(ncidF, 'Tyz', TyzFid)
   end if
   ierr=nf90_get_var(ncidS,TxxSid,ExxS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyySid,EyyS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TzzSid,EzzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxySid,ExyS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TxzSid,ExzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidS,TyzSid,EyzS,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxxFid,ExxR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyyFid,EyyR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TzzFid,EzzR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxyFid,ExyR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TxzFid,ExzR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
   ierr=nf90_get_var(ncidF,TyzFid,EyzR,(/bsubs,1/),(/bsubc,knt/),(/bsubt,1/))
end if

! load media
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt) 
! convert to stain
  call stress2strain(ExxS,EyyS,EzzS,ExyS,ExzS,EyzS,lambda,mu,bsubc)
  call stress2strain(ExxR,EyyR,EzzR,ExyR,ExzR,EyzR,lambda,mu,bsubc)

! cal kernel
   call cal_kernel(V,U,                &
        ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
        ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
        lambda,mu,                     &
        nt1,nt2,kdt,bsubc,             &
        Kap,Kaq,Kbp,Kbq)
! put
   call nfseis_put(ncid,kapid,Kap,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,kaqid,Kaq,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,kbpid,Kbp,bsubs,bsubc,bsubt)
   call nfseis_put(ncid,kbqid,Kbq,bsubs,bsubc,bsubt)

end do
end do
end do

if (snap_tcnt(kid)>=knt) then
   call nfseis_close(ncidS)
   call nfseis_close(ncidF)
end if

   call nfseis_close(ncid)

end if kernode_if

!-----------------------------------------------------------------------------
call media_destroy
call dealloc_local

#ifdef USEMPI
call MPI_BARRIER(SWMPI_COMM,ierr)
call MPI_FINALIZE(ierr)
#endif

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine read_tomo_conf(fnm_conf)
character (len=*),intent(in) :: fnm_conf
character (len=SEIS_STRLEN) :: fnm_filter
integer fid,n

fid=5002
open(fid,file=trim(fnm_conf),status='old')
call string_conf(fid,1,'MAIN_CONF',2,fnm_main_conf)

call string_conf(fid,1,'record_file',2,fnm_record)
call string_conf(fid,1,'time_window',2,kwin(1))
call string_conf(fid,1,'time_window',3,kwin(2))
call string_conf(fid,1,'filter_file',2,fnm_filter)

call string_conf(fid,1,'WAVE_ROOT',2,pnm_wave)
call string_conf(fid,1,'SGT_ROOT',2,pnm_sgt)
call string_conf(fid,1,"sgt_source_M0",2,src_m0)

call string_conf(fid,1,'KERNEL_ROOT',2,pnm_ker)
call string_conf(fid,1,'kernel_prefix',2,ker_prefix)
call string_conf(fid,1,'kernel_snap_id',2,kid)

do n=1,SEIS_GEO
 call string_conf(fid,1,'process_unit',n+1,blklen(n))
end do
close(fid)

open(fid,file=trim(fnm_record),status='old')
  read(fid,*) knt
  call alloc_seismo(knt)
  do n=1,knt
     read(fid,*) T(n),V(n),U(n)
  end do
close(fid)

flag_filter=.false.
if (trim(fnm_filter)/='none') then
   flag_filter=.true.
   open(fid,file=trim(fnm_filter),status='old')
   read(fid,*) NOD1
   NOD=NOD1-1; NEL=knt+2*NOD1+NOD
   allocate(filt_b(NOD1)); filt_b=0.0;
   allocate(filt_a(NOD1)); filt_a=0.0
   allocate(vecx(NEL)); vecx=0.0
   allocate(vecy(NEL)); vecy=0.0
   allocate(vecw(NEL)); vecw=0.0
   do n=1,NOD1
      read(fid,*) filt_a(n),filt_b(n)
   end do
   close(fid)
end if

end subroutine read_tomo_conf

subroutine alloc_seismo(nt)
  integer,intent(in) :: nt
  allocate(V(nt)); V=0.0
  allocate(U(nt)); U=0.0
  allocate(T(nt)); T=0.0
end subroutine alloc_seismo

subroutine ker_para_reinit
  where (blklen==-1)
     blklen=snap_subc(:,kid)  
  end where
  blknum=(snap_subc(:,kid)+blklen-1)/blklen
end subroutine ker_para_reinit

subroutine alloc_local(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0
  allocate(ExxR(ki,kj,kk,kt)); ExxR=0.0
  allocate(EyyR(ki,kj,kk,kt)); EyyR=0.0
  allocate(EzzR(ki,kj,kk,kt)); EzzR=0.0
  allocate(ExyR(ki,kj,kk,kt)); ExyR=0.0
  allocate(ExzR(ki,kj,kk,kt)); ExzR=0.0
  allocate(EyzR(ki,kj,kk,kt)); EyzR=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0

  allocate(Ka(kt),Kb(kt));Ka=0.0;Kb=0.0
  allocate(E11S(kt)); E11S=0.0
  allocate(E22S(kt)); E22S=0.0
  allocate(E33S(kt)); E33S=0.0
  allocate(E12S(kt)); E12S=0.0
  allocate(E13S(kt)); E13S=0.0
  allocate(E23S(kt)); E23S=0.0
  allocate(E11R(kt)); E11R=0.0
  allocate(E22R(kt)); E22R=0.0
  allocate(E33R(kt)); E33R=0.0
  allocate(E12R(kt)); E12R=0.0
  allocate(E13R(kt)); E13R=0.0
  allocate(E23R(kt)); E23R=0.0
end subroutine alloc_local

subroutine alloc_media_local(ki,kj,kk)
  integer,intent(in) :: ki,kj,kk
  allocate(lambda(ki,kj,kk)); lambda=0.0
  allocate(mu(ki,kj,kk)); mu=0.0
end subroutine alloc_media_local

subroutine dealloc_local
  if (allocated(V)) deallocate(V)
  if (allocated(U)) deallocate(U)
  if (allocated(T)) deallocate(T)
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
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)
end subroutine dealloc_local

subroutine realloc_local(ki,kj,kk,kt)
  integer,intent(in) :: ki,kj,kk,kt
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
  if (allocated(Kap)) deallocate(Kap)
  if (allocated(Kaq)) deallocate(Kaq)
  if (allocated(Kbp)) deallocate(Kbp)
  if (allocated(Kbq)) deallocate(Kbq)

  allocate(ExxS(ki,kj,kk,kt)); ExxS=0.0
  allocate(EyyS(ki,kj,kk,kt)); EyyS=0.0
  allocate(EzzS(ki,kj,kk,kt)); EzzS=0.0
  allocate(ExyS(ki,kj,kk,kt)); ExyS=0.0
  allocate(ExzS(ki,kj,kk,kt)); ExzS=0.0
  allocate(EyzS(ki,kj,kk,kt)); EyzS=0.0
  allocate(ExxR(ki,kj,kk,kt)); ExxR=0.0
  allocate(EyyR(ki,kj,kk,kt)); EyyR=0.0
  allocate(EzzR(ki,kj,kk,kt)); EzzR=0.0
  allocate(ExyR(ki,kj,kk,kt)); ExyR=0.0
  allocate(ExzR(ki,kj,kk,kt)); ExzR=0.0
  allocate(EyzR(ki,kj,kk,kt)); EyzR=0.0

  allocate(Kap(ki,kj,kk   )); Kap=0.0 
  allocate(Kaq(ki,kj,kk   )); Kaq=0.0 
  allocate(Kbp(ki,kj,kk   )); Kbp=0.0 
  allocate(Kbq(ki,kj,kk   )); Kbq=0.0 
end subroutine realloc_local

subroutine cal_kernel(V,U,                &
           ExxR,EyyR,EzzR,ExyR,ExzR,EyzR, &
           ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
           lam,miu,                       &
           nt1,nt2,dt,subc,Kap,Kaq,Kbp,Kbq)

integer,intent(in) :: nt1,nt2
real(SP),dimension(:),intent(in) :: V,U
real(SP),dimension(:,:,:,:),intent(in) ::         &
     ExxS,EyyS,EzzS,ExyS,ExzS,EyzS, &
     ExxR,EyyR,EzzR,ExyR,ExzR,EyzR
real(SP),dimension(:,:,:),intent(in) :: miu,lam
real(SP),dimension(:,:,:),intent(out) :: Kap,Kaq,Kbp,Kbq
integer,dimension(SEIS_GEO),intent(in) :: subc
real(SP),intent(in) :: dt

integer :: m,n,i,j,k,ts,te
real(SP) :: Vprho,Vsrho
real(kind=SP) :: x1,x2

if (flag_filter) then
   ts=1; te=knt
else
   ts=nt1; te=nt2
end if

do k=1,subc(3)
do j=1,subc(2)
do i=1,subc(1)
   E11R(1:te)=ExxR(i,j,k,1:te);E22R=EyyR(i,j,k,1:te);E33R=EzzR(i,j,k,1:te)
   E12R(1:te)=ExyR(i,j,k,1:te);E13R=ExzR(i,j,k,1:te);E23R=EyzR(i,j,k,1:te)
   E11S(1:te)=ExxS(i,j,k,1:te);E22S=EyyS(i,j,k,1:te);E33S=EzzS(i,j,k,1:te)
   E12S(1:te)=ExyS(i,j,k,1:te);E13S=ExzS(i,j,k,1:te);E23S=EyzS(i,j,k,1:te)

do n=ts,te
   x1=0.0_SP; x2=0.0_SP
do m=1,n-1
   x1=x1+1.0_SP*dt*(E11R(m)+E22R(m)+E33R(m))*(E11S(n-m)+E22S(n-m)+E33S(n-m))
   x2=x2+1.0_SP*dt*( E11R(m)*E11S(n-m)+E22R(m)*E22S(n-m)+E33R(m)*E33S(n-m)     &
               +2.0_SP*E12R(m)*E12S(n-m)+2.0_SP*E13R(m)*E13S(n-m)+2.0_SP*E23R(m)*E23S(n-m) )
end do
   Ka(n)=x1; Kb(n)=2.0_SP*(x2-x1)
end do

if (flag_filter) then
   vecx(NOD+1:NOD+te)=Ka
   call filtfilt(NOD,NEL,filt_b,filt_a,vecx,vecy,vecw)
   Ka(nt1:nt2)=vecy(NOD+nt1:NOD+nt2)
   
   vecx(NOD+1:NOD+te)=Kb
   call filtfilt(NOD,NEL,filt_b,filt_a,vecx,vecy,vecw)
   Kb(nt1:nt2)=vecy(NOD+nt1:NOD+nt2)
end if

! Forming kernels for tau_p and tau_q.
Kap(i,j,k)=0.5*dt*(V(nt1)*Ka(nt1)+V(nt2)*Ka(nt2))
Kaq(i,j,k)=0.5*dt*(U(nt1)*Ka(nt1)+U(nt2)*Ka(nt2))
Kbp(i,j,k)=0.5*dt*(V(nt1)*Kb(nt1)+V(nt2)*Kb(nt2))
Kbq(i,j,k)=0.5*dt*(U(nt1)*Kb(nt1)+U(nt2)*Kb(nt2))
do n=nt1+1,nt2-1
   Kap(i,j,k)=Kap(i,j,k)+V(n)*Ka(n)*dt
   Kaq(i,j,k)=Kaq(i,j,k)+U(n)*Ka(n)*dt
   Kbp(i,j,k)=Kbp(i,j,k)+V(n)*Kb(n)*dt
   Kbq(i,j,k)=Kbq(i,j,k)+U(n)*Kb(n)*dt
end do

Vprho= lam(i,j,k)+2.0*miu(i,j,k); Vsrho= miu(i,j,k)
Kap(i,j,k)= Kap(i,j,k)*2.0*Vprho
Kaq(i,j,k)=-Kaq(i,j,k)*2.0*Vprho
Kbp(i,j,k)= Kbp(i,j,k)*2.0*Vsrho
Kbq(i,j,k)=-Kbq(i,j,k)*2.0*Vsrho

end do
end do
end do
end subroutine cal_kernel

subroutine stress2strain(Exx,Eyy,Ezz,Exy,Exz,Eyz,lam,miu,scl)
real(SP),dimension(:,:,:,:),intent(inout) :: Exx,Eyy,Ezz,Exy,Exz,Eyz
real(SP),dimension(:,:,:),intent(in) :: lam,miu
integer,dimension(SEIS_GEO),intent(in) :: scl
real(SP) E1,E2,E3,E0
integer i,j,k,n,n0

n0=size(Exx,4)
do n=1,n0
do k=1,scl(3)
do j=1,scl(2)
do i=1,scl(1)
   E1=(lam(i,j,k)+miu(i,j,k))/(miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   E2=-lam(i,j,k)/(2.0*miu(i,j,k)*(3.0*lam(i,j,k)+2.0*miu(i,j,k)))
   E3=1.0/miu(i,j,k)
   E0=E2*(Exx(i,j,k,n)+Eyy(i,j,k,n)+Ezz(i,j,k,n))
   Exx(i,j,k,n)=E0-(E2-E1)*Exx(i,j,k,n)
   Eyy(i,j,k,n)=E0-(E2-E1)*Eyy(i,j,k,n)
   Ezz(i,j,k,n)=E0-(E2-E1)*Ezz(i,j,k,n)
   Exy(i,j,k,n)=0.5*E3*Exy(i,j,k,n)
   Exz(i,j,k,n)=0.5*E3*Exz(i,j,k,n)
   Eyz(i,j,k,n)=0.5*E3*Eyz(i,j,k,n)
end do
end do
end do
end do
end subroutine stress2strain

subroutine filtfilt(NOD,NEL,b,a,x,y,w)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real,dimension(:) :: x,y,w

x(1:NOD)=0.0; x(NEL-2*(NOD+1)+1:NEL)=0.0
call filter(NOD,NEL,b,a,x,y)
w=y(NEL:1:-1)
call filter(NOD,NEL,b,a,w,y)
y=y(NEL:1:-1)
end subroutine filtfilt

subroutine filter(NOD,NEL,b,a,x,y)
integer,intent(in) :: NOD,NEL
real(DP),dimension(:),intent(in) :: b,a
real,dimension(:) :: x,y
integer :: n,m

y(1:NOD)=0.0
do n=NOD+1,NEL
   y(n)=b(1)*x(n)
   do m=2,NOD+1
      y(n)=y(n)+b(m)*x(n-m+1)-a(m)*y(n-m+1)
   end do
   y(n)=y(n)/a(1)
end do
end subroutine filter

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
#ifdef USEMPI
  integer :: ierr
#endif
  print *, trim(msg)
#ifdef USEMPI
  call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
  stop 1
#endif
end subroutine error_except

end program tomo_kernel

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:

