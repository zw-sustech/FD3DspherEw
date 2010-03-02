module media_mod

! This module contains the variables for 3D medium
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG

!*****************************************************************************
!
! $Date$
! $Revision$
! $LastChangedBy$
!
!*****************************************************************************

use constants_mod
use string_mod
use math_mod
use para_mod
use mpi_mod
use nfseis_mod

implicit none
private
public ::         &
  media_fnm_init, &
  media_fnm_get,  &
  media_destroy,  &
  media_alloc,    &
  media_import

!-----------------------------------------------------------------------------

real(SP),dimension(:,:,:),allocatable,public :: rho,buoy
real(SP),dimension(:,:,:),allocatable,public :: C11,C12,C13,C14,C15,C16
real(SP),dimension(:,:,:),allocatable,public ::     C22,C23,C24,C25,C26
real(SP),dimension(:,:,:),allocatable,public ::         C33,C34,C35,C36
real(SP),dimension(:,:,:),allocatable,public ::             C44,C45,C46
real(SP),dimension(:,:,:),allocatable,public ::                 C55,C56
real(SP),dimension(:,:,:),allocatable,public ::                     C66
real(SP),dimension(:,:,:),allocatable,public :: Qs
real(SP),public :: QsF0,QsINF
character (len=SEIS_STRLEN),public ::       &
     fnm_media_conf, pnm_media
integer :: ierr

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

!*************************************************************************
!*                            alloc and dealloc                          *
!*************************************************************************
subroutine media_alloc
  allocate( rho(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  rho=0.0_SP

  allocate( C13(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  C13=0.0_SP
  allocate( C66(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  C66=0.0_SP
#if defined AnisVTI || defined AnisGene
  allocate( C11(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  C11=0.0_SP
  allocate( C33(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  C33=0.0_SP
  allocate( C44(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  C44=0.0_SP
#endif
#ifdef AnisGene
 !allocate( C11(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C11=0.0_SP
  allocate( C12(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C12=0.0_SP
 !allocate( C13(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C13=0.0_SP
  allocate( C14(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C14=0.0_SP
  allocate( C15(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C15=0.0_SP
  allocate( C16(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C16=0.0_SP
  allocate( C22(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C22=0.0_SP
  allocate( C23(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C23=0.0_SP
  allocate( C24(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C24=0.0_SP
  allocate( C25(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C25=0.0_SP
  allocate( C26(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C26=0.0_SP
 !allocate( C33(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C33=0.0_SP
  allocate( C34(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C34=0.0_SP
  allocate( C35(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C35=0.0_SP
  allocate( C36(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C36=0.0_SP
 !allocate( C44(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C44=0.0_SP
  allocate( C45(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C45=0.0_SP
  allocate( C46(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C46=0.0_SP
  allocate( C55(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C55=0.0_SP
  allocate( C56(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C56=0.0_SP
 !allocate( C66(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); C66=0.0_SP
#endif
#ifdef WITHQS
  allocate( Qs(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Qs=0.0_SP
#endif
end subroutine media_alloc
subroutine media_destroy
  if (allocated(rho)) deallocate( rho )
  if (allocated(Qs )) deallocate( Qs )
  if (allocated(C11)) deallocate(C11)
  if (allocated(C12)) deallocate(C12)
  if (allocated(C13)) deallocate(C13)
  if (allocated(C14)) deallocate(C14)
  if (allocated(C15)) deallocate(C15)
  if (allocated(C16)) deallocate(C16)
  if (allocated(C22)) deallocate(C22)
  if (allocated(C23)) deallocate(C23)
  if (allocated(C24)) deallocate(C24)
  if (allocated(C25)) deallocate(C25)
  if (allocated(C26)) deallocate(C26)
  if (allocated(C33)) deallocate(C33)
  if (allocated(C34)) deallocate(C34)
  if (allocated(C35)) deallocate(C35)
  if (allocated(C36)) deallocate(C36)
  if (allocated(C44)) deallocate(C44)
  if (allocated(C45)) deallocate(C45)
  if (allocated(C46)) deallocate(C46)
  if (allocated(C55)) deallocate(C55)
  if (allocated(C56)) deallocate(C56)
  if (allocated(C66)) deallocate(C66)
end subroutine media_destroy

!*************************************************************************
!*                             media io                                  *
!*************************************************************************
subroutine media_fnm_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
    call string_conf(fid,1,'MEDIA_CONF',2,fnm_media_conf)
    call string_conf(fid,1,'MEDIA_ROOT',2,pnm_media)
  close(fid)
end subroutine media_fnm_init

function media_fnm_get(n_i,n_j,n_k) result(filenm)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_media)//'/'//'media'//'_'//set_mpi_subfix(n_i,n_j,n_k)//'.nc'
end function media_fnm_get

subroutine media_import(n_i,n_j,n_k)
  integer,intent(in) :: n_i,n_j,n_k
  character (len=SEIS_STRLEN) :: filenm
  integer,dimension(SEIS_GEO) :: subs,subc,subt
  subs=(/ 1,1,1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
  filenm=media_fnm_get(n_i,n_j,n_k)
  call nfseis_varget( filenm, 'rho', rho, subs,subc,subt)
  call nfseis_varget( filenm, 'C66', C66, subs,subc,subt)

  call nfseis_varget( filenm, 'C13', C13, subs,subc,subt)
#if defined AnisVTI || defined AnisGene
  call nfseis_varget( filenm, 'C11', C11, subs,subc,subt)
  call nfseis_varget( filenm, 'C33', C33, subs,subc,subt)
  call nfseis_varget( filenm, 'C44', C44, subs,subc,subt)
#endif
#ifdef AnisGene
 !call nfseis_varget( filenm, 'C11', C11, subs,subc,subt)
  call nfseis_varget( filenm, 'C12', C12, subs,subc,subt)
 !call nfseis_varget( filenm, 'C13', C13, subs,subc,subt)
  call nfseis_varget( filenm, 'C14', C14, subs,subc,subt)
  call nfseis_varget( filenm, 'C15', C15, subs,subc,subt)
  call nfseis_varget( filenm, 'C16', C16, subs,subc,subt)
  call nfseis_varget( filenm, 'C22', C22, subs,subc,subt)
  call nfseis_varget( filenm, 'C23', C23, subs,subc,subt)
  call nfseis_varget( filenm, 'C24', C24, subs,subc,subt)
  call nfseis_varget( filenm, 'C25', C25, subs,subc,subt)
  call nfseis_varget( filenm, 'C26', C26, subs,subc,subt)
 !call nfseis_varget( filenm, 'C33', C33, subs,subc,subt)
  call nfseis_varget( filenm, 'C34', C34, subs,subc,subt)
  call nfseis_varget( filenm, 'C35', C35, subs,subc,subt)
  call nfseis_varget( filenm, 'C36', C36, subs,subc,subt)
 !call nfseis_varget( filenm, 'C44', C44, subs,subc,subt)
  call nfseis_varget( filenm, 'C45', C45, subs,subc,subt)
  call nfseis_varget( filenm, 'C46', C46, subs,subc,subt)
  call nfseis_varget( filenm, 'C55', C55, subs,subc,subt)
  call nfseis_varget( filenm, 'C56', C56, subs,subc,subt)
 !call nfseis_varget( filenm, 'C66', C66, subs,subc,subt)
#endif
#ifdef WITHQS
  call nfseis_varget( filenm, 'Qs', Qs, subs,subc,subt)
  call nfseis_attget( filenm, 'QsF0', QsF0)
  call nfseis_attget( filenm, 'QsINF', QsINF)
#endif
end subroutine media_import
  
end module media_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
