program seis3d_media

!-------------------------------------------------------------------------------
! Description:
!-------------------------------------------------------------------------------
!
! This program inits the medium of the 3D geology model
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG
!
!-------------------------------------------------------------------------------
! Time stamp, log, version:
!-------------------------------------------------------------------------------
!
! $Date$
! $Revision$
! $LastChangedBy$
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Include file:
!-------------------------------------------------------------------------------
#include "mod_macdrp.h"

!-------------------------------------------------------------------------------
! Dependent modules:
!-------------------------------------------------------------------------------
    use constants_mod
    use string_mod
    use math_mod
    use para_mod
    use mpi_mod
    use nfseis_mod
    use grid_mod
    use media_mod
    use custom_mod
    use io_mod
#ifdef MediaMPI
    use mpi
#endif

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
    implicit none
    
    integer,parameter ::           &
         SIG_ANIS_ISO       =100,  &
         SIG_ANIS_VTI       =200,  &
         SIG_ANIS_TTI       =300,  &
         SIG_ANIS_TRICLINIC =310,  &
         SIG_STIF_VEL       =  1,  &
         SIG_STIF_TEN       =  2,  &
         SIG_ISO_TEN = SIG_ANIS_ISO + SIG_STIF_TEN, &
         SIG_ISO_VEL = SIG_ANIS_ISO + SIG_STIF_VEL, &
         SIG_VTI_TEN = SIG_ANIS_VTI + SIG_STIF_TEN, &
         SIG_VTI_VEL = SIG_ANIS_VTI + SIG_STIF_VEL, &
         SIG_TTI_TEN = SIG_ANIS_TTI + SIG_STIF_TEN, &
         SIG_TTI_VEL = SIG_ANIS_TTI + SIG_STIF_VEL
    
    integer :: NGRIDX,NGRIDY,NGRIDZ
    
    type STRUCT_1D
         logical :: yes
         integer :: nk
         character (len=SEIS_STRLEN) :: fnm
         character (len=SEIS_STRLEN) :: filetype
         real(SP),dimension(:),pointer :: z
         real(SP),dimension(:),pointer :: d
         real(SP),dimension(:),pointer :: Vp,Vs,Dp
         real(SP),dimension(:),pointer :: Vph,Vsh,eta
         real(SP),dimension(:),pointer :: Qp,Qs
         real(SP) :: QsF0,QsINF
    end type STRUCT_1D
    
    type STRUCT_INTERFACE
         logical :: yes
         character (len=SEIS_STRLEN) :: fnm
         character (len=SEIS_STRLEN) :: filetype
         integer :: ni,nj,nk
         real(SP),dimension(:),pointer :: x,y
         real(SP),dimension(:,:,:),pointer :: z
         real(SP),dimension(:,:,:),pointer :: h
         real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
         real(SP),dimension(:,:),pointer :: Vph,Vsh,eta
         real(SP),dimension(:,:),pointer :: theta,phi
         real(SP),dimension(:,:,:),pointer :: C
         real(SP),dimension(:,:),pointer :: Qs
         real(SP) :: QsF0,QsINF
    end type STRUCT_INTERFACE
    
    type STRUCT_LAYERED
         logical :: yes
         character (len=SEIS_STRLEN) :: fnm
         character (len=SEIS_STRLEN) :: filetype
         integer :: ni,nj,nk
         real(SP),dimension(:),pointer :: x,y
         real(SP),dimension(:,:,:),pointer :: d
         real(SP),dimension(:,:,:),pointer :: h
         real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
         real(SP),dimension(:,:),pointer :: Vph,Vsh,eta
         real(SP),dimension(:,:),pointer :: theta,phi
         real(SP),dimension(:,:,:),pointer :: C
         real(SP),dimension(:,:),pointer :: Qs
         real(SP) :: QsF0,QsINF
    end type STRUCT_LAYERED
    
    type STRUCT_COMPOSITE
         logical :: yes
         character (len=SEIS_STRLEN) :: fnm
         character (len=SEIS_STRLEN) :: filetype
         integer :: ni,nj,nk
         real(SP),dimension(:),pointer :: x,y
         real(SP),dimension(:,:,:),pointer :: h
         real(SP),dimension(:,:,:,:),pointer :: Vp,Vs,Dp
         real(SP),dimension(:,:,:,:),pointer :: Vph,Vsh,eta
         real(SP),dimension(:,:,:,:),pointer :: theta,phi
         real(SP),dimension(:,:,:,:,:),pointer :: C
         real(SP),dimension(:,:,:,:),pointer :: Qs
         real(SP),dimension(:),pointer :: Vp_poly_d,Vs_poly_d,Dp_poly_d,Qs_poly_d
         real(SP),dimension(:),pointer :: Vph_poly_d,Vsh_poly_d,eta_poly_d
         real(SP),dimension(:),pointer :: theta_poly_d,phi_poly_d
         real(SP),dimension(:,:),pointer :: C_poly_d
         real(SP) :: QsF0,QsINF
    end type STRUCT_COMPOSITE
    
!-- for structure sampled on volume nodes
!-- since the large data size, only read a sub-domain containing this thread
    type STRUCT_VOLUME
         logical :: yes
         character (len=SEIS_STRLEN) :: fnm
         character (len=SEIS_STRLEN) :: filetype
         integer :: imax,jmax,kmax
         integer :: ni,nj,nk
         real(SP) :: rz0
         real(SP),dimension(:),pointer :: gx,gy,grz
         real(SP),dimension(:),pointer :: x,y,rz
         real(SP),dimension(:,:,:),pointer :: Vp,Vs,Dp
         real(SP),dimension(:,:,:),pointer :: Vph,Vsh,eta
         real(SP),dimension(:,:,:),pointer :: theta,phi
         real(SP),dimension(:,:,:,:),pointer :: C
         real(SP),dimension(:,:,:),pointer :: Qs
         real(SP) :: QsF0,QsINF
    end type STRUCT_VOLUME
    
!-- for structure sampled on horizontal parallelpipe nodes, and vertically
!-- interfaces, with polynomial coefficients and power
!-- since the large data size, only read a sub-domain containing this thread
    type STRUCT_VERPOLY
         logical :: yes
         character (len=SEIS_STRLEN) :: fnm
         character (len=SEIS_STRLEN) :: filetype
         integer :: imax,jmax,kmax,nmax
         integer :: ni,nj,nk
         real(SP) :: rz0
         real(SP),dimension(:),pointer :: gx,gy
         real(SP),dimension(:),pointer :: x,y
         real(SP),dimension(:,:,:),pointer :: rz
         real(SP),pointer :: Vp_poly_c(:,:,:,:),Vp_poly_d(:)
         real(SP),pointer :: Vs_poly_c(:,:,:,:),Vs_poly_d(:)
         real(SP),pointer :: Dp_poly_c(:,:,:,:),Dp_poly_d(:)
         real(SP),pointer :: Qs_poly_c(:,:,:,:),Qs_poly_d(:)
         real(SP) :: QsF0,QsINF
    end type STRUCT_VERPOLY
    
    type (STRUCT_1D) :: L1D
    type (STRUCT_INTERFACE) :: LA
    type (STRUCT_LAYERED) :: L3D
    type (STRUCT_COMPOSITE) :: LC
    type (STRUCT_VOLUME) :: BV,PV
    type (STRUCT_VERPOLY) :: BP,PP
    integer :: flag_anis_type,flag_stif_type
    
    real(SP) :: ztopo
    
    real(DP),dimension(:,:,:),allocatable :: gx,gy,gz
    real(SP) :: dtmax,dtmaxVp,dtmaxL
    integer,dimension(SEIS_GEO) :: dtindx,dtnode
    
    integer :: n_i,n_j,n_k
#ifdef MediaMPI
    integer :: ierr,ROOT
    real(SP),dimension(2) :: dtin,dtout
#endif

!-------------------------------------------------------------------------------
! Program Entries:
!-------------------------------------------------------------------------------

#ifdef MediaMPI
    call MPI_INIT(ierr)
#endif

    call get_conf_name(fnm_conf)

    call swmpi_init(fnm_conf)
    call para_init(fnm_conf)

#ifdef MediaMPI
    call swmpi_cart_creat
    call swmpi_reinit_para
    call swmpi_datatype
#else
    call swmpi_set_gindx(0,0,0)
#endif

    call grid_fnm_init(fnm_conf)
    call media_fnm_init(fnm_conf)
    call media_alloc
    call grid_alloc
    
    call alloc_local

!-------------------------------------------------------------------------------
!-- init media configure and read input media --
!-------------------------------------------------------------------------------
#ifdef MediaMPI
    if (masternode)  &
#endif
       print *, 'init media ...'
    call init_media(fnm_media_conf)

!-------------------------------------------------------------------------------
!-- perform medium volumeric integration --
!-------------------------------------------------------------------------------
    dtmax=1.0e10

#ifdef MediaMPI
    if (masternode)  &
#endif
       print *, 'calculate effective media ...'

#ifdef MediaMPI
    n_i=thisid(1); n_j=thisid(2); n_k=thisid(3)
#else
    do n_i=0,dims(1)-1
    do n_j=0,dims(2)-1
    do n_k=dims(3)-1,0,-1
       write(*,"(i10,2(i3),a,3(i3))") n_i,n_j,n_k, ' of ',dims
       call swmpi_change_fnm(n_i,n_j,n_k)
       call swmpi_set_gindx(n_i,n_j,n_k)
#endif

       !=== load axis from the free surface block to get ztopo ==
       call grid_coord_import(n_i,n_j,dims(3)-1); ztopo=z(nk2)
       
       !=== load grid axis for this block ===
       call grid_coord_import(n_i,n_j,n_k)
    
       !=== read background volume data in this block ===
       if (BV%yes) call volume_read(BV)
       !=== read perturbed volume data in this block ===
       if (PV%yes) call volume_read(PV)
    
       !=== read background verpoly data in this block ===
       if (BP%yes) call verpoly_read(BP)
       !=== read perturbed verpoly data in this block ===
       if (PP%yes) call verpoly_read(PP)
    
       !=== calculate effective medium parameters for this block ===
       call effmedia_eval
    
       !=== extend medium parameters to ghost grid points ===
       call media_extend

#ifdef MediaMPI
       if (masternode) write(*,*) "exchange media on boundary stencil ..."
       call media_exchange
#endif

#ifdef MediaMPI
       if (masternode)  &
#endif
          print *, '  export media...'
       !=== create and output ===
       call media_create(n_i,n_j,n_k)
    
       !=== calculate maximum time step based on distance estimation ===
       call media_stept_check(n_i,n_j,n_k)
    
#ifndef MediaMPI
    end do
    end do
    end do

   !=== change info between blocks ===
   print *, "exchange media on boundary stencil ..."
   call media_exchange
#endif

!-------------------------------------------------------------------------------
!-- gather largest stept among blocks if mpi is used, then print --
!-------------------------------------------------------------------------------
#ifdef MediaMPI
    dtin=(/ dtmax,real(myid,SP) /)
    ROOT=0
#ifdef DataTypeDouble
    call MPI_REDUCE(dtin,dtout,1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,ROOT,SWMPI_COMM,ierr)
#else
    call MPI_REDUCE(dtin,dtout,1,MPI_2REAL,MPI_MINLOC,ROOT,SWMPI_COMM,ierr)
#endif
    dtmax=dtout(1)
    if (myid==ROOT) then
    print *, "Maximum allowed time step is", dtmax, ' in thread', nint(dtout(2))
    if (dtmax<stept) then
       print *, "Serious Error: stept>dtmax", stept,dtmax
       print *, " occurs on thread", nint(dtout(2))
       !call MPI_ABORT(SWMPI_COMM,1,ierr)
    end if
    end if
#else
    print *, "Maximum allowed time step is", dtmax
    write(*,"(a,3i5,a,3i5)") "located on", dtindx,' in thread', dtnode
    print *, " Vp and dL are:", dtmaxVp,dtmaxL
    if (dtmax<stept) then
       print *, "Serious Error: stept>dtmax", stept,dtmax
       !stop 1
    end if
#endif

!-------------------------------------------------------------------------------
!-- code finished, clean workspace
!-------------------------------------------------------------------------------
    call media_destroy
    call grid_dealloc
    call local_destroy

#ifdef MediaMPI
    call MPI_BARRIER(SWMPI_COMM,ierr)
    call MPI_FINALIZE(ierr)
#endif

!===============================================================================
contains
!===============================================================================

!===============================================================================
subroutine init_media(fnm_conf)
!===============================================================================
    character (len=*),intent(in) :: fnm_conf
    character (len=SEIS_STRLEN) :: str
    integer :: fid
    
    fid=1001
    !== opon the conf file ==
    open(fid,file=trim(fnm_conf),status="old")
    
    !== read volume-integration point number ==
    call string_conf(fid,1,'effective_sampling_point',2,NGRIDX)
    call string_conf(fid,1,'effective_sampling_point',3,NGRIDY)
    call string_conf(fid,1,'effective_sampling_point',4,NGRIDZ)
    
    !== initialize all instances ==
    L1D%yes=.false.; L3D%yes=.false.;
    LA%yes=.false.; LC%yes=.false.;
    BP%yes=.false.; BV%yes=.false.;
    PP%yes=.false.; PV%yes=.false.;
    
    !== read background model information ==
    call string_conf(fid,1,'background_type',2,str)
    select case (trim(str))
    case ('cart1d')
         L1D%yes=.true.
         call string_conf(fid,1,'background_format',2,L1D%filetype)
         call string_conf(fid,1,'background_filename',2,L1D%fnm)
         call layer1d_read(L1D)
    case ('interface')
         LA%yes=.true.
         call string_conf(fid,1,'background_format',2,LA%filetype)
         call string_conf(fid,1,'background_filename',2,LA%fnm)
         call interface_read(LA)
    case ('layered')
         L3D%yes=.true.
         call string_conf(fid,1,'background_format',2,L3D%filetype)
         call string_conf(fid,1,'background_filename',2,L3D%fnm)
         call layered_read(L3D)
    case ('composite')
         LC%yes=.true.
         call string_conf(fid,1,'background_format',2,LC%filetype)
         call string_conf(fid,1,'background_filename',2,LC%fnm)
         call composite_read(LC)
    case ('volume')
         BV%yes=.true.
         call string_conf(fid,1,'background_format',2,BV%filetype)
         call string_conf(fid,1,'background_filename',2,BV%fnm)
         call volume_init(BV)
    case ('verpoly')
         BP%yes=.true.
         call string_conf(fid,1,'background_format',2,BP%filetype)
         call string_conf(fid,1,'background_filename',2,BP%fnm)
         call verpoly_init(BP)
    end select
    
    !== read perturbed model information ==
    call string_conf(fid,1,'perturbed_type',2,str)
    select case (trim(str))
    case ('volume')
         PV%yes=.true.
         call string_conf(fid,1,'perturbed_format',2,PV%filetype)
         call string_conf(fid,1,'perturbed_filename',2,PV%fnm)
         call volume_init(PV)
    case ('verpoly')
         PP%yes=.true.
         call string_conf(fid,1,'perturbed_format',2,PP%filetype)
         call string_conf(fid,1,'perturbed_filename',2,PP%fnm)
         call verpoly_init(PP)
    case ('none')
    end select
    
    !== close the conf file ==
    close(fid)
    
    !if (W%yes .and. L%yes) then
    !   print *, "cann't deal with poly and layered media at same time currently"
    !   stop 1
    !end if
!===============================================================================
end subroutine init_media
!===============================================================================

!-------------------------------------------------------------------------------
subroutine local_destroy
!-------------------------------------------------------------------------------
    call layer1d_destroy
    call layered_destroy
    call interface_destroy
    call composite_destroy
    call volume_destroy
    call verpoly_destroy
!-------------------------------------------------------------------------------
end subroutine local_destroy
!-------------------------------------------------------------------------------

!===============================================================================
!-------------------------------  1d structure --------------------------------
!===============================================================================

!-------------------------------------------------------------------------------
subroutine layer1d_read(L)
!-------------------------------------------------------------------------------
!== 1D structure is modified from normal mode code input file.
    type(STRUCT_1D) :: L
    integer :: fid,kmax,k
    integer :: ifanis, ifdeck, nic, noc
    real(SP) :: corrfreq
    real(DP) :: d0
    character (len=SEIS_STRLEN) :: str
    
    if (trim(L%filetype)/='ascii') then
       call error_except('only ascii input accepted for cart1d type')
    end if
    
    fid=2001
    open(fid,file=trim(L%fnm),status="old")
    
    read(fid,"(a)") str
    read(fid,*) ifanis, corrfreq, ifdeck
    read(fid,*) kmax, nic, noc
    
    if (ifanis==1) then
       flag_anis_type=SIG_ANIS_VTI
    else
       flag_anis_type=SIG_ANIS_ISO
    end if
    flag_stif_type=stif_name2id('velocity')
    
    L%nk=kmax
    
    allocate(L%z(kmax))
    allocate(L%d(kmax))
    allocate(L%Vp(kmax))
    allocate(L%Vs(kmax))
    allocate(L%Dp(kmax))
    allocate(L%Qp(kmax))
    allocate(L%Qs(kmax))
    allocate(L%Vph(kmax))
    allocate(L%Vsh(kmax))
    allocate(L%eta(kmax))
    
    do k=1,kmax
       read(fid,*) L%z(k),L%Dp(k),L%Vp(k),L%Vs(k),L%Qp(k),L%Qs(k), &
                   L%Vph(k),L%Vsh(k),L%eta(k)
    end do
    
    close(fid)
    
    ! check
    if ( any(L%Dp<=0) .or. any(L%Vp<=0) .or. any(L%Vs<0) ) then
       print *, 'media parameter is negtive'
       print *, 'Vp=',L%Vp
       print *, 'Vs=',L%Vs
       print *, 'rho=',L%Dp
       call error_except('cart type read failed')
    end if
    if ( any(L%Vp**2<=2.0*L%Vs**2) ) then
       print *, 'Vp**2 < 2Vs**2'
       print *, 'Vp=',L%Vp
       print *, 'Vs=',L%Vs
       call error_except('cart type read failed')
    end if
    
    d0=max(L%z(1),L%z(kmax))
    do k=1,kmax
       L%d(k)=d0-L%z(k)
    end do
    
    if (L%d(1)>L%d(kmax)) then
       L%d  =L%d  (kmax:1:-1)
       L%Vp =L%Vp (kmax:1:-1)
       L%Vs =L%Vs (kmax:1:-1)
       L%Dp =L%Dp (kmax:1:-1)
       L%Vph=L%Vph(kmax:1:-1)
       L%Vsh=L%Vsh(kmax:1:-1)
       L%eta=L%eta(kmax:1:-1)
    end if

#ifdef WITHQS
    QsF0=L%QsF0; QsINF=L%QsINF
#endif
!-------------------------------------------------------------------------------
end subroutine layer1d_read
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine layer1d_retrieve(d0,Vp,Vs,Dp,Vph,Vsh,eta)
!-------------------------------------------------------------------------------
!-- Purpose: get parameters values from 1D structure for given d0. 
!-- 1D structure is modified from normal mode code input file.
!-- In 1D structure, the parameter vector is order from surface to depth.
!-- Two entries for an interface, same depth but different medium parameters.

    real(SP),intent(in) :: d0   != input: depth of position
    real(SP),dimension(:),intent(out) :: &
        Vp,Vs,Dp,Vph,Vsh,eta    != output: medium parameters
    
    real(SP) :: d,h,L1,L2
    integer :: k1,k2,k0
    
    !== init index surrounding d0 to the begin and end of the vector ==
    k1=1; k2=L1D%nk
    
    !== if d0 is at the surface ==
    if (abs(L1D%d(k1)-d0)<SEIS_ZERO .or. d0<L1D%d(k1)) then
       Vp (:)=L1D%Vp (k1); Vs (:)=L1D%Vs (k1); Dp (:)=L1D%Dp (k1)
       Vph(:)=L1D%Vph(k1); Vsh(:)=L1D%Vsh(k1); eta(:)=L1D%eta(k1)
    !== if d0 at the bottom ==
    elseif (abs(L1D%d(k2)-d0)<SEIS_ZERO) then
       !-- d(k2) is an interface --
       if (abs(L1D%d(k2-1)-L1D%d(k2))<SEIS_ZERO) then
          Vp (1)=L1D%Vp (k2-1);Vp (2)=L1D%Vp (k2);
          Vs (1)=L1D%Vs (k2-1);Vs (2)=L1D%Vs (k2)
          Dp (1)=L1D%Dp (k2-1);Dp (2)=L1D%Dp (k2)
          Vph(1)=L1D%Vph(k2-1);Vph(2)=L1D%Vph(k2);
          Vsh(1)=L1D%Vsh(k2-1);Vsh(2)=L1D%Vsh(k2)
          eta(1)=L1D%eta(k2-1);eta(2)=L1D%eta(k2)
       !-- otherwise --
       else
          Vp (:)=L1D%Vp (k2); Vs (:)=L1D%Vs (k2); Dp (:)=L1D%Dp (k2)
          Vph(:)=L1D%Vph(k2); Vsh(:)=L1D%Vsh(k2); eta(:)=L1D%eta(k2)
       end if
    !== if d0 below the lowest depth in the 1D vector ==
    elseif (L1D%d(k2)<d0) then
       Vp (:)=L1D%Vp (k2); Vs (:)=L1D%Vs (k2); Dp (:)=L1D%Dp (k2)
       Vph(:)=L1D%Vph(k2); Vsh(:)=L1D%Vsh(k2); eta(:)=L1D%eta(k2)
    !== otherwise narrow k1 k2 to the nearest index ==
    else
        !-- search k1 k2 by 2-divided method --
        do while (k2-k1>1)
           k0=(k2+k1)/2
           !-- exit if d0 exactly locates at half of k1 and k2 --
           if (abs(L1D%d(k0)-d0)<SEIS_ZERO) then
              k1=k0;k2=k0; exit
           !-- if d0 above d(k0), lift k2 to k0 --
           elseif (L1D%d(k0)>d0) then
              k2=k0;
           !-- if d0 below d(k), put down k1 to k0 --
           else
              k1=k0
           end if
        end do
    
        !-- if d0 exactly locates at some index --
        if (k1==k2) then
           if (abs(L1D%d(k1-1)-L1D%d(k1))<SEIS_ZERO) then
              Vp (1)=L1D%Vp (k1-1); Vp (2)=L1D%Vp (k1);
              Vs (1)=L1D%Vs (k1-1); Vs (2)=L1D%Vs (k1)
              Dp (1)=L1D%Dp (k1-1); Dp (2)=L1D%Dp (k1)
              Vph(1)=L1D%Vph(k1-1); Vph(2)=L1D%Vph(k1);
              Vsh(1)=L1D%Vsh(k1-1); Vsh(2)=L1D%Vsh(k1)
              eta(1)=L1D%eta(k1-1); eta(2)=L1D%eta(k1)
           elseif (abs(L1D%d(k1+1)-L1D%d(k1))<SEIS_ZERO) then
              Vp (1)=L1D%Vp (k1); Vp (2)=L1D%Vp (k1+1);
              Vs (1)=L1D%Vs (k1); Vs (2)=L1D%Vs (k1+1)
              Dp (1)=L1D%Dp (k1); Dp (2)=L1D%Dp (k1+1)
              Vph(1)=L1D%Vph(k1); Vph(2)=L1D%Vph(k1+1);
              Vsh(1)=L1D%Vsh(k1); Vsh(2)=L1D%Vsh(k1+1)
              eta(1)=L1D%eta(k1); eta(2)=L1D%eta(k1+1)
           else
              Vp (:)=L1D%Vp (k1); Vs (:)=L1D%Vs (k1); Dp (:)=L1D%Dp (k1)
              Vph(:)=L1D%Vph(k1); Vsh(:)=L1D%Vsh(k1); eta(:)=L1D%eta(k1)
           end if
        !-- interpolate if d0 is between k1 and k2
        else
           d=L1D%d(k2); !-d0;
           h=L1D%d(k2)-L1D%d(k1)
           L1=(d-d0)/h; L2=1.0-L1
           Vp (:)=L1D%Vp (k1)*L1+L1D%Vp (k2)*L2
           Dp (:)=L1D%Dp (k1)*L1+L1D%Dp (k2)*L2
           Vph(:)=L1D%Vph(k1)*L1+L1D%Vph(k2)*L2

           !Vs (:)=L1D%Vs (k1)*L1+L1D%Vs (k2)*L2
           !Vsh(:)=L1D%Vsh(k1)*L1+L1D%Vsh(k2)*L2
           !eta(:)=L1D%eta(k1)*L1+L1D%eta(k2)*L2

           if (L1D%Vs(k1) < SEIS_ZERO .or. L1D%Vs(k2) < SEIS_ZERO) then
            Vs (1)=L1D%Vs (k1); Vs (2)=L1D%Vs (k2); 
           else
            Vs (:)=L1D%Vs (k1)*L1+L1D%Vs (k2)*L2
           end if

           if (L1D%Vsh(k1) < SEIS_ZERO .or. L1D%Vsh(k2) < SEIS_ZERO) then
            Vsh (1)=L1D%Vsh (k1); Vsh (2)=L1D%Vsh (k2); 
           else
            Vsh (:)=L1D%Vsh (k1)*L1+L1D%Vsh (k2)*L2
           end if

           if (L1D%eta(k1) < SEIS_ZERO .or. L1D%eta(k2) < SEIS_ZERO) then
            eta (1)=L1D%eta (k1); eta (2)=L1D%eta (k2); 
           else
            eta (:)=L1D%eta (k1)*L1+L1D%eta (k2)*L2
           end if
        end if
    end if
!-------------------------------------------------------------------------------
end subroutine layer1d_retrieve
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine layer1d_destroy
!-------------------------------------------------------------------------------
    if (associated(L1D%z  )) deallocate(L1D%z  )
    if (associated(L1D%d  )) deallocate(L1D%d  )
    if (associated(L1D%Dp )) deallocate(L1D%Dp )
    if (associated(L1D%Vp )) deallocate(L1D%Vp )
    if (associated(L1D%Vs )) deallocate(L1D%Vs )
    if (associated(L1D%Qp )) deallocate(L1D%Qp )
    if (associated(L1D%Qs )) deallocate(L1D%Qs )
    if (associated(L1D%Vph)) deallocate(L1D%Vph)
    if (associated(L1D%Vsh)) deallocate(L1D%Vsh)
    if (associated(L1D%eta)) deallocate(L1D%eta)
!-------------------------------------------------------------------------------
end subroutine layer1d_destroy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine interface_read(L)
!-------------------------------------------------------------------------------
!-- Purpose: Read interface structure input file.
!-- Interface is a contact surface between two media.
!-- The model is given by the interface axis and the medium parameters 
!-- above and below it along with a vertical linear variation constrain.
!-- The interface is ordered from the free surface to depth.
!-- The free surface is taken as an interface between air and the first 
!-- layer. The medium below the last interface is an half space.

    type(STRUCT_INTERFACE) :: L
    integer :: fid,imax,jmax,kmax,i,j,k,n
    character (len=SEIS_STRLEN) :: str
    
!-- only ascii file valid --
    if (trim(L%filetype)/='ascii') then
       call error_except('only ascii input accepted for interface type')
    end if
    
!-- open the interface file to read --
    fid=2001
    open(fid,file=trim(L%fnm),status="old")
    
!-- anisotrop_type --
    call string_conf(fid,1,'anisotropy_type',2,str)
    flag_anis_type=anis_name2id(trim(str))
    
!-- stiffness_type --
    call string_conf(fid,1,'stiffness_type',2,str)
    flag_stif_type=stif_name2id(trim(str))
    
!-- read number of interface into kmax --
    call string_conf(fid,1,'number_of_interface',2,kmax)
    if (kmax<1) then
       call error_except('at least 1 interface should exist for interface type')
    end if
    
!-- number of lateral sampling point --
    call string_conf(fid,1,'horizontal_sampling',2,imax)
    call string_conf(fid,1,'horizontal_sampling',3,jmax)
    
!-- keep discrete sampling information --
    L%ni=imax; L%nj=jmax; L%nk=kmax
    
!-------- read in interface --------
    allocate(L%x(imax))
    allocate(L%y(jmax))
    allocate(L%z(imax,jmax,kmax))
    allocate(L%h(imax,jmax,kmax-1))
    
    call string_conf(fid,1,'<anchor_interface>',1,str)
    do j=1,jmax
    do i=1,imax     
        read(fid,*) L%x(i),L%y(j),( L%z(i,j,k), k=1,kmax )
    end do 
    end do
    
!-- convert interface to layer thickness --
    do k=1,kmax-1
       L%h(:,:,k)=L%z(:,:,k)-L%z(:,:,k+1)
    end do
    
!-- check layer thinckness non-negative --
    do k=1,kmax-1
       if (any(L%h(:,:,k)<0.0)) then
          print *, 'thickness of layer',k,' should be not less than 0'
          print *, minloc(L%h(:,:,k))
          call error_except("thickness btween two interfaces shouldn' be less than 0")
       end if
    end do
    
!-- convert degree to radian --
    L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP

!-------- read in medium parameters --------

#ifdef WITHQS
!-- Grave's approximate anelastic parameters --
    call string_conf(fid,1,'QsF0',2,L%QsF0)
    call string_conf(fid,1,'QsINF',2,L%QsINF)
    allocate(L%Qs(2,kmax))
    QsF0=L%QsF0; QsINF=L%QsINF
    call error_except("WITHQS is ongoing")
#endif

    !-- allocate rho --
    allocate(L%Dp(2,kmax)) !-- 2 is for above and below --
    
    !-- if input stiffness tensor --
    if (flag_stif_type==SIG_STIF_TEN) then
       allocate(L%C(21,2,kmax)); L%C=0.0_SP
       allocate(L%theta(2,kmax)); L%theta=0.0_SP
       allocate(L%phi(2,kmax)); L%phi=0.0_SP
    
       !-- ISO
       if (flag_anis_type==SIG_ANIS_ISO) then
          call string_conf(fid,1,'<anchor_ISO_tensor>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), L%C(3,1,k), L%C(21,1,k)
             read(fid,*) L%Dp(2,k), L%C(3,2,k), L%C(21,2,k)
          end do
       end if
       !-- VTI
       if (flag_anis_type==SIG_ANIS_VTI) then
          call string_conf(fid,1,'<anchor_VTI_tensor>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), L%C(1,1,k),L%C(3,1,k),L%C(12,1,k),L%C(16,1,k),L%C(21,1,k)
             read(fid,*) L%Dp(2,k), L%C(1,2,k),L%C(3,2,k),L%C(12,2,k),L%C(16,2,k),L%C(21,2,k)
          end do
       end if
       !-- TTI
       if (flag_anis_type==SIG_ANIS_TTI) then
          call string_conf(fid,1,'<anchor_TTI_tensor>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), L%C(1,1,k),L%C(3,1,k),L%C(12,1,k),L%C(16,1,k),L%C(21,1,k),L%theta(1,k),L%phi(1,k)
             read(fid,*) L%Dp(2,k), L%C(1,2,k),L%C(3,2,k),L%C(12,2,k),L%C(16,2,k),L%C(21,2,k),L%theta(2,k),L%phi(2,k)
          end do
          L%theta=L%theta*PI/180.0_DP
          L%phi=L%phi*PI/180.0_DP
       end if
       !-- TRICLINIC
       if (flag_anis_type==SIG_ANIS_TRICLINIC) then
          call string_conf(fid,1,'<anchor_TRICLINIC_tensor>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), ( L%C(n,1,k), n=1,21 )
             read(fid,*) L%Dp(2,k), ( L%C(n,2,k), n=1,21 )
          end do
       end if
    end if
  
!-- velocity parameters
    if (flag_stif_type==SIG_STIF_VEL) then
       allocate(L%Vp (2,kmax)); L%Vp =0.0_SP
       allocate(L%Vs (2,kmax)); L%Vs =0.0_SP
       allocate(L%Vph(2,kmax)); L%Vph=0.0_SP
       allocate(L%Vsh(2,kmax)); L%Vsh=0.0_SP
       allocate(L%eta(2,kmax)); L%eta=0.0_SP
       allocate(L%theta(2,kmax)); L%theta=0.0_SP
       allocate(L%phi(2,kmax)); L%phi=0.0_SP
    
       !-- ISO
       if (flag_anis_type==SIG_ANIS_ISO) then
          call string_conf(fid,1,'<anchor_ISO_velocity>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), L%Vp(1,k), L%Vs(1,k)
             read(fid,*) L%Dp(2,k), L%Vp(2,k), L%Vs(2,k)
          end do
       end if
       !-- VTI
       if (flag_anis_type==SIG_ANIS_VTI) then
          call string_conf(fid,1,'<anchor_VTI_velocity>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), L%Vp(1,k),L%Vs(1,k),L%Vph(1,k),L%Vph(1,k),L%eta(1,k)
             read(fid,*) L%Dp(2,k), L%Vp(2,k),L%Vs(2,k),L%Vph(2,k),L%Vph(2,k),L%eta(2,k)
          end do
       end if
       !-- TTI
       if (flag_anis_type==SIG_ANIS_TTI) then
          call string_conf(fid,1,'<anchor_TTI_tensor>',1,str)
          do k=1,kmax
             read(fid,*) L%Dp(1,k), L%Vp(1,k),L%Vs(1,k),L%Vph(1,k),L%Vph(1,k),L%eta(1,k),L%phi(1,k),L%theta(1,k)
             read(fid,*) L%Dp(2,k), L%Vp(2,k),L%Vs(2,k),L%Vph(2,k),L%Vph(2,k),L%eta(2,k),L%phi(2,k),L%theta(2,k)
          end do
          L%theta=L%theta*PI/180.0_DP
          L%phi=L%phi*PI/180.0_DP
       end if
       !-- TRICLINIC
       if (flag_anis_type==SIG_ANIS_TRICLINIC) then
          call error_except("no velocity equivalent for TRICLINIC medium")
       end if
    end if

    ! check
    !do k=1,kmax
    !do i=1,2
    !   if (k==1 .and. i==1) then
    !      if (L%Dp(i,k)<0 .or. L%Vp(i,k)<0 .or. L%Vs(i,k)<0 ) then
    !          call error_except('model parameter should >=0 for (1,1)')
    !      end if
    !   elseif ( L%Dp(i,k)<=0 .or. L%Vp(i,k)<=0 .or. L%Vs(i,k)<0 ) then
    !      print *, 'media parameter is negtive'
    !      print *, 'Vp=',L%Vp
    !      print *, 'Vs=',L%Vs
    !      print *, 'Dp=',L%Dp
    !      call error_except('interface type read failed')
    !   elseif ( L%Vp(i,k)**2<=2.0*L%Vs(i,k)**2 ) then
    !      print *, 'Vp**2 < 2Vs**2'
    !      print *, 'Vp=',L%Vp
    !      print *, 'Vs=',L%Vs
    !      call error_except('interface type read failed')
    !   end if
    !end do
    !end do
    
    !-- close the interface file --
    close(fid)
!-------------------------------------------------------------------------------
end subroutine interface_read
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine interface_destroy
!-------------------------------------------------------------------------------
    if (associated(LA%x  )  ) deallocate(LA%x  )
    if (associated(LA%y  )  ) deallocate(LA%y  )
    if (associated(LA%z  )  ) deallocate(LA%z  )
    if (associated(LA%h  )  ) deallocate(LA%h  )
    if (associated(LA%Dp )  ) deallocate(LA%Dp )
    if (associated(LA%Vp )  ) deallocate(LA%Vp )
    if (associated(LA%Vs )  ) deallocate(LA%Vs )
    !if (associated(LA%Qp )  ) deallocate(LA%Qp )
    if (associated(LA%Qs )  ) deallocate(LA%Qs )
    if (associated(LA%Vph)  ) deallocate(LA%Vph)
    if (associated(LA%Vsh)  ) deallocate(LA%Vsh)
    if (associated(LA%eta)  ) deallocate(LA%eta)
    if (associated(LA%theta)) deallocate(LA%theta)
    if (associated(LA%phi)  ) deallocate(LA%phi)
    if (associated(LA%C)    ) deallocate(LA%C)
!-------------------------------------------------------------------------------
end subroutine interface_destroy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine interface_retrieve(x0,y0,z0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: x0,y0,z0,dz2
    real(SP),dimension(:),intent(out) :: Vp,Vs,Dp,Vph,Vsh,eta,theta,phi
    real(SP),dimension(:,:),intent(out) :: CMat
    
    real(SP) :: x1,y1,z1,z2
    real(SP) :: d,h,L1,L2
    
    integer :: n,m
    integer :: i1,i2,j1,j2,k1,k2,k0,n1,n2
    
    call indx_locate_1d(x0,LA%x,i1,i2,x1)
    call indx_locate_1d(y0,LA%y,j1,j2,y1)
    
    do n=LA%nk,1,-1
    
       z1=interp_2d(LA%x(i1:i2),LA%y(j1:j2),LA%z(i1:i2,j1:j2,n),2,2,x1,y1)
       if (n==LA%nk) then
          h=1.0e3
       else
          h =interp_2d(LA%x(i1:i2),LA%y(j1:j2),LA%h(i1:i2,j1:j2,n),2,2,x1,y1)
       end if
    
       if (h<=SEIS_ZERO .and. n>1) cycle
    
       if (h<=SEIS_ZERO .and. n==1) then
          if (LA%Vp(1,n)>SEIS_ZERO) then
             n1=1;k1=1;n2=1;k2=1
          end if
          ! get value
          Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
             Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(1)=LA%Vph(n1,k1); Vph(2)=LA%Vph(n2,k2)
                Vsh(1)=LA%Vsh(n1,k1); Vsh(2)=LA%Vsh(n2,k2)
                eta(1)=LA%eta(n1,k1); eta(2)=LA%eta(n2,k2)
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(1)=LA%theta(n1,k1); theta(2)=LA%theta(n2,k2)
                phi(1)=LA%phi(n1,k1); phi(2)=LA%phi(n2,k2)
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(1,m)=LA%C(m,n1,k1); CMat(2,m)=LA%C(m,n2,k2)
             end do
          end if
          exit
       end if
    
       if (abs(z1-z0)<=dz2/5.0) then  !  .and. h>SEIS_ZERO
          k2=n; n2=2
          if (n==1) then
             if (LA%Vp(1,n)>SEIS_ZERO) then
                k1=n; n1=1
             else
                k1=n; n1=2
             end if
          end if
          do m=n-1,1,-1
             h=interp_2d(LA%x(i1:i2),LA%y(j1:j2),LA%h(i1:i2,j1:j2,m),2,2,x1,y1)
             if (h>SEIS_ZERO) then
                 k1=m+1; n1=1
                 exit
             elseif (m==1) then
                if (LA%Vp(1,m)>SEIS_ZERO) then
                   k1=m; n1=1
                else
                   k1=n; n1=2
                end if
             end if
          end do
          ! get value
          Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
             Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(1)=LA%Vph(n1,k1); Vph(2)=LA%Vph(n2,k2)
                Vsh(1)=LA%Vsh(n1,k1); Vsh(2)=LA%Vsh(n2,k2)
                eta(1)=LA%eta(n1,k1); eta(2)=LA%eta(n2,k2)
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(1)=LA%theta(n1,k1); theta(2)=LA%theta(n2,k2)
                phi(1)=LA%phi(n1,k1); phi(2)=LA%phi(n2,k2)
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(1,m)=LA%C(m,n1,k1); CMat(2,m)=LA%C(m,n2,k2)
             end do
          end if
          exit
       end if
    
       if (z0<z1) then
          if (n==LA%nk) then
             Dp(:)=LA%Dp(2,n)
             if (flag_stif_type==SIG_STIF_VEL) then
                Vp(:)=LA%Vp(2,n)
                Vs(:)=LA%Vs(2,n)
                if (flag_anis_type>=SIG_ANIS_VTI) then
                   Vph(:)=LA%Vph(2,n)
                   Vsh(:)=LA%Vsh(2,n)
                   eta(:)=LA%eta(2,n)
                end if
                if (flag_anis_type>=SIG_ANIS_TTI) then
                   theta(:)=LA%theta(2,n)
                   phi(:)=LA%phi(2,n)
                end if
             end if
             if (flag_stif_type==SIG_STIF_TEN) then
                do m=1,21
                   Cmat(:,m)=LA%C(m,2,n)
                end do
             end if
          else
             L2=(z1-z0)/h; L1=1.0-L2
             Dp(:)=LA%Dp(2,n)*L1+LA%Dp(1,n+1)*L2
             if (flag_stif_type==SIG_STIF_VEL) then
                Vp(:)=LA%Vp(2,n)*L1+LA%Vp(1,n+1)*L2
                Vs(:)=LA%Vs(2,n)*L1+LA%Vs(1,n+1)*L2
                if (flag_anis_type>=SIG_ANIS_VTI) then
                   Vph(:)=LA%Vph(2,n)*L1+LA%Vph(1,n+1)*L2
                   Vsh(:)=LA%Vsh(2,n)*L1+LA%Vsh(1,n+1)*L2
                   eta(:)=LA%eta(2,n)*L1+LA%eta(1,n+1)*L2
                end if
                if (flag_anis_type>=SIG_ANIS_TTI) then
                   theta(:)=LA%theta(2,n)*L1+LA%theta(1,n+1)*L2
                   phi(:)=LA%phi(2,n)*L1+LA%phi(1,n+1)*L2
                end if
             end if
             if (flag_stif_type==SIG_STIF_TEN) then
                do m=1,21
                   CMat(:,m)=LA%C(m,2,n)*L1+LA%C(m,1,n+1)*L2
                end do
             end if
          end if
          exit
       end if
    
       k1=n; k2=n; n1=2; n2=2
    
       if (n==1) then
          if (LA%Vp(1,n)>SEIS_ZERO) then
             n1=1;k1=1;n2=1;k2=1
          end if
          Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
             Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(1)=LA%Vph(n1,k1); Vph(2)=LA%Vph(n2,k2)
                Vsh(1)=LA%Vsh(n1,k1); Vsh(2)=LA%Vsh(n2,k2)
                eta(1)=LA%eta(n1,k1); eta(2)=LA%eta(n2,k2)
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(1)=LA%theta(n1,k1); theta(2)=LA%theta(n2,k2)
                phi(1)=LA%phi(n1,k1); phi(2)=LA%phi(n2,k2)
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(1,m)=LA%C(m,n1,k1); CMat(2,m)=LA%C(m,n2,k2)
             end do
          end if
          exit
       end if
    
    end do
!-------------------------------------------------------------------------------
end subroutine interface_retrieve
!-------------------------------------------------------------------------------

!===============================================================================
! -----------------------------  layered structure ----------------------------
!===============================================================================

!-------------------------------------------------------------------------------
subroutine layered_read(L)
!-------------------------------------------------------------------------------
    type(STRUCT_LAYERED) :: L
    integer :: fid,imax,jmax,kmax,i,j,k,n
    character (len=SEIS_STRLEN) :: str,layer_type
    
    if (trim(L%filetype)/='ascii') then
       call error_except('only ascii input accepted for layered type')
    end if
    
    fid=2001
    open(fid,file=trim(L%fnm),status="old")
    call string_conf(fid,1,'number_of_layer',2,kmax)
    if (kmax<1) then
       call error_except('at least 1 layer should exist for layered type')
    end if
    
#ifdef WITHQS
    call string_conf(fid,1,'QsF0',2,L%QsF0)
    call string_conf(fid,1,'QsINF',2,L%QsINF)
#endif
    
    !call string_conf(fid,1,'layer_sealevel',2,layer_sealevel)
    call string_conf(fid,1,'horizontal_sampling',2,imax)
    call string_conf(fid,1,'horizontal_sampling',3,jmax)
    
    call string_conf(fid,1,'anisotropy_type',2,str)
    flag_anis_type=anis_name2id(trim(str))
    
    call string_conf(fid,1,'stiffness_type',2,str)
    flag_stif_type=stif_name2id(trim(str))
    
    L%ni=imax; L%nj=jmax; L%nk=kmax
    
!-- read in interface --
    allocate(L%x(imax))
    allocate(L%y(jmax))
    allocate(L%d(imax,jmax,kmax))
    allocate(L%h(imax,jmax,kmax))
    
    call string_conf(fid,1,'layer_meaning',2,layer_type)
    call string_conf(fid,1,'<anchor_layer>',1,str)
    select case (trim(layer_type))
    case ('depth')
      do j=1,jmax
      do i=1,imax     
          read(fid,*) L%x(i),L%y(j),( L%d(i,j,k), k=1,kmax )
      end do 
      end do
      L%h(:,:,1)=L%d(:,:,1)
      do k=2,kmax
         L%h(:,:,k)=L%d(:,:,k)-L%d(:,:,k-1)
      end do
    case ('thickness')
      do j=1,jmax
      do i=1,imax     
          read(fid,*) L%x(i),L%y(j),( L%h(i,j,k), k=1,kmax )
      end do 
      end do
      L%d(:,:,1)=L%h(:,:,1)
      do k=2,kmax
         L%d(:,:,k)=L%d(:,:,k-1)+L%h(:,:,k)
      end do
    case default
         call error_except("layer_meaning can only take depth or thickness")
    end select
    
    do k=1,kmax-1
       if (any(L%h(:,:,k)<0.0)) then
          print *, minloc(L%h(:,:,k))
          call error_except('thickness of layer should not be less than 0')
       end if
    end do
    if (any(L%h(:,:,kmax)<SEIS_ZERO)) then
       print *, minloc(L%h(:,:,kmax))
       call error_except('thickness of lowest layer should larger than 0')
    end if
    
    L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP
    
!-- read in medium parameters --
    !-- rho
    allocate(L%Dp(2,kmax))
    call string_conf(fid,1,'<anchor_rho>',1,str)
    do k=1,kmax
       read(fid,*) L%Dp(:,k)
    end do
    
#ifdef WITHQS
    !-- Qs
    allocate(L%Qs(2,kmax))
    call string_conf(fid,1,'<anchor_Qs>',1,str)
    do k=1,kmax
       read(fid,*) L%Qs(:,k)
    end do
#endif
    
!-- stiffness tensor
    if (flag_stif_type==SIG_STIF_TEN) then
       allocate(L%C(21,2,kmax)); L%C=0.0_SP
       allocate(L%theta(2,kmax)); L%theta=0.0_SP
       allocate(L%phi(2,kmax)); L%phi=0.0_SP
    
       !-- ISO
       if (flag_anis_type==SIG_ANIS_ISO) then
          call string_conf(fid,1,'<anchor_stiffness_ISO>',1,str)
          do k=1,kmax
             read(fid,*) L%C(3,1,k), L%C(21,1,k)
             read(fid,*) L%C(3,2,k), L%C(21,2,k)
          end do
       end if
       !-- VTI
       if (flag_anis_type==SIG_ANIS_VTI) then
          call string_conf(fid,1,'<anchor_stiffness_VTI>',1,str)
          do k=1,kmax
             read(fid,*) L%C(1,1,k),L%C(3,1,k),L%C(12,1,k),L%C(16,1,k),L%C(21,1,k)
             read(fid,*) L%C(1,2,k),L%C(3,2,k),L%C(12,2,k),L%C(16,2,k),L%C(21,2,k)
          end do
       end if
       !-- TTI
       if (flag_anis_type==SIG_ANIS_TTI) then
          call string_conf(fid,1,'<anchor_stiffness_TTI>',1,str)
          do k=1,kmax
             read(fid,*) L%C(1,1,k),L%C(3,1,k),L%C(12,1,k),L%C(16,1,k),L%C(21,1,k),L%theta(1,k),L%phi(1,k)
             read(fid,*) L%C(1,2,k),L%C(3,2,k),L%C(12,2,k),L%C(16,2,k),L%C(21,2,k),L%theta(2,k),L%phi(2,k)
          end do
          L%theta=L%theta*PI/180.0_DP
          L%phi=L%phi*PI/180.0_DP
       end if
       !-- TRICLINIC
       if (flag_anis_type==SIG_ANIS_TRICLINIC) then
          call string_conf(fid,1,'<anchor_stiffness_TRICLINIC>',1,str)
          do k=1,kmax
             read(fid,*) ( L%C(n,1,k), n=1,21 )
             read(fid,*) ( L%C(n,2,k), n=1,21 )
          end do
       end if
    end if
    
!-- velocity parameters
    if (flag_stif_type==SIG_STIF_VEL) then
       allocate(L%Vp (2,kmax)); L%Vp =0.0_SP
       allocate(L%Vs (2,kmax)); L%Vs =0.0_SP
       allocate(L%Vph(2,kmax)); L%Vph=0.0_SP
       allocate(L%Vsh(2,kmax)); L%Vsh=0.0_SP
       allocate(L%eta(2,kmax)); L%eta=0.0_SP
       allocate(L%theta(2,kmax)); L%theta=0.0_SP
       allocate(L%phi(2,kmax)); L%phi=0.0_SP
    
       !-- ISO
       if (flag_anis_type>=SIG_ANIS_ISO) then
          call string_conf(fid,1,'<anchor_Vp>',1,str)
          do k=1,kmax
             read(fid,*) L%Vp(1,k), L%Vp(2,k)
          end do
          call string_conf(fid,1,'<anchor_Vs>',1,str)
          do k=1,kmax
             read(fid,*) L%Vs(1,k), L%Vs(2,k)
          end do
       end if
       !-- VTI
       if (flag_anis_type>=SIG_ANIS_VTI) then
          call string_conf(fid,1,'<anchor_Vph>',1,str)
          do k=1,kmax
             read(fid,*) L%Vph(1,k), L%Vph(2,k)
          end do
          call string_conf(fid,1,'<anchor_Vsh>',1,str)
          do k=1,kmax
             read(fid,*) L%Vsh(1,k), L%Vsh(2,k)
          end do
          call string_conf(fid,1,'<anchor_eta>',1,str)
          do k=1,kmax
             read(fid,*) L%eta(1,k), L%eta(2,k)
          end do
       end if
       !-- TTI
       if (flag_anis_type>=SIG_ANIS_TTI) then
          call string_conf(fid,1,'<anchor_theta>',1,str)
          do k=1,kmax
             read(fid,*) L%theta(1,k), L%theta(2,k)
          end do
          call string_conf(fid,1,'<anchor_phi>',1,str)
          do k=1,kmax
             read(fid,*) L%phi(1,k), L%phi(2,k)
          end do
          L%theta=L%theta*PI/180.0_DP
          L%phi=L%phi*PI/180.0_DP
       end if
       !-- TRICLINIC
       if (flag_anis_type==SIG_ANIS_TRICLINIC) then
          call error_except("no velocity equivalent for TRICLINIC medium")
       end if
    end if
    
    !! check
    !if ( any(L%Dp<=0) .or. any(L%Vp<=0) .or. any(L%Vs<0) ) then
    !   print *, 'media parameter is negtive'
    !   print *, 'Vp=',L%Vp
    !   print *, 'Vs=',L%Vs
    !   print *, 'Dp=',L%Dp
    !   call error_except('layered type read failed')
    !end if
    !if ( any(L%Vp**2<=2.0*L%Vs**2) ) then
    !   print *, 'Vp**2 < 2Vs**2'
    !   print *, 'Vp=',L%Vp
    !   print *, 'Vs=',L%Vs
    !   call error_except('layered type read failed')
    !end if
    
    close(fid)
    
#ifdef WITHQS
    QsF0=L%QsF0; QsINF=L%QsINF
#endif
!-------------------------------------------------------------------------------
end subroutine layered_read
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine layered_destroy
!-------------------------------------------------------------------------------
    if (associated(L3D%x  )  ) deallocate(L3D%x  )
    if (associated(L3D%y  )  ) deallocate(L3D%y  )
    if (associated(L3D%d  )  ) deallocate(L3D%d  )
    if (associated(L3D%h  )  ) deallocate(L3D%h  )
    if (associated(L3D%Dp )  ) deallocate(L3D%Dp )
    if (associated(L3D%Vp )  ) deallocate(L3D%Vp )
    if (associated(L3D%Vs )  ) deallocate(L3D%Vs )
    !if (associated(L3D%Qp )  ) deallocate(L3D%Qp )
    if (associated(L3D%Qs )  ) deallocate(L3D%Qs )
    if (associated(L3D%Vph)  ) deallocate(L3D%Vph)
    if (associated(L3D%Vsh)  ) deallocate(L3D%Vsh)
    if (associated(L3D%eta)  ) deallocate(L3D%eta)
    if (associated(L3D%theta)) deallocate(L3D%theta)
    if (associated(L3D%phi)  ) deallocate(L3D%phi)
    if (associated(L3D%C)    ) deallocate(L3D%C)
!-------------------------------------------------------------------------------
end subroutine layered_destroy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine layered_retrieve(x0,y0,d0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: x0,y0,d0,dz2
    real(SP),dimension(:),intent(out) :: Vp,Vs,Dp,Vph,Vsh,eta,theta,phi
    real(SP),dimension(:,:),intent(out) :: CMat
    
    real(SP) :: x1,y1,z1,z2
    real(SP) :: d,h,L1,L2
    
    integer :: n,m
    integer :: i1,i2,j1,j2,k1,k2,k0,n1,n2
    
    d=0.0
    call indx_locate_1d(x0,L3D%x,i1,i2,x1)
    call indx_locate_1d(y0,L3D%y,j1,j2,y1)
    do n=1,L3D%nk
       h=interp_2d(L3D%x(i1:i2),L3D%y(j1:j2),L3D%h(i1:i2,j1:j2,n),2,2,x1,y1)
    
       if (h<=SEIS_ZERO) cycle
    
       d=d+h
       if (abs(d-d0)<=dz2/5.0) then
          if (n==L3D%nk) call error_except('please increase the thickness of the lowest layer')
          n1=n;
          do m=n+1,L3D%nk
             h=interp_2d(L3D%x(i1:i2),L3D%y(j1:j2),L3D%h(i1:i2,j1:j2,m),2,2,x1,y1)
             if (h>SEIS_ZERO) then
                 n2=m
                 exit
             end if
          end do
          ! get value
          Dp(1)=L3D%Dp(2,n1); Dp(2)=L3D%Dp(1,n2)
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(1)=L3D%Vp(2,n1); Vp(2)=L3D%Vp(1,n2)
             Vs(1)=L3D%Vs(2,n1); Vs(2)=L3D%Vs(1,n2)
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(1)=L3D%Vph(2,n1); Vph(2)=L3D%Vph(1,n2)
                Vsh(1)=L3D%Vsh(2,n1); Vsh(2)=L3D%Vsh(1,n2)
                eta(1)=L3D%eta(2,n1); eta(2)=L3D%eta(1,n2)
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(1)=L3D%theta(2,n1); theta(2)=L3D%theta(1,n2)
                phi(1)=L3D%phi(2,n1); phi(2)=L3D%phi(1,n2)
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(1,m)=L3D%C(m,2,n1); CMat(2,m)=L3D%C(m,1,n2)
             end do
          end if
          exit
       elseif (d>d0) then
          L1=(d-d0)/h; L2=1.0-L1
          ! get value
          Dp(:)=L3D%Dp(1,n)*L1+L3D%Dp(2,n)*L2
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(:)=L3D%Vp(1,n)*L1+L3D%Vp(2,n)*L2
             Vs(:)=L3D%Vs(1,n)*L1+L3D%Vs(2,n)*L2
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(:)=L3D%Vph(1,n)*L1+L3D%Vph(2,n)*L2
                Vsh(:)=L3D%Vsh(1,n)*L1+L3D%Vsh(2,n)*L2
                eta(:)=L3D%eta(1,n)*L1+L3D%eta(2,n)*L2
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(:)=L3D%theta(1,n)*L1+L3D%theta(2,n)*L2
                phi(:)=L3D%phi(1,n)*L1+L3D%phi(2,n)*L2
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(:,m)=L3D%C(m,1,n)*L1+L3D%C(m,2,n)*L2
             end do
          end if
          exit
       elseif (n==L3D%nk) then
          call error_except('please increase the thickness of the lowest layer')
       end if
    end do
!-------------------------------------------------------------------------------
end subroutine layered_retrieve
!-------------------------------------------------------------------------------

!===============================================================================
! -----------------------------  composite structure ----------------------------
!===============================================================================

!-------------------------------------------------------------------------------
subroutine composite_read(L)
!-------------------------------------------------------------------------------
    type(STRUCT_COMPOSITE) :: L
    integer :: imax,jmax,kmax
    character (len=SEIS_STRLEN) :: ctype
    
    if (trim(L%filetype)/='nc') then
       call error_except('only nc input accepted for composite type')
    end if
    
    call nfseis_diminfo(L%fnm,'theta',imax)
    call nfseis_diminfo(L%fnm,'phi',jmax)
    call nfseis_diminfo(L%fnm,'layer',kmax)
    
    call nfseis_attget(L%fnm,'anisotrop_type',ctype)
    flag_anis_type=anis_name2id(trim(adjustl(ctype)))
    
    call nfseis_attget(L%fnm,'stiffness_type',ctype)
    flag_stif_type=stif_name2id(trim(adjustl(ctype)))
    
    L%ni=imax; L%nj=jmax; L%nk=kmax
    
!-- coordinate
    allocate(L%x(imax))
    allocate(L%y(jmax))
    allocate(L%h(imax,jmax,kmax))
    call nfseis_varget(L%fnm,'theta',L%x,(/1/),(/imax/),(/1/))
    call nfseis_varget(L%fnm,'phi',L%y,(/1/),(/jmax/),(/1/))
    call nfseis_varget(L%fnm,'thickness',L%h,(/1,1,1/),(/imax,jmax,kmax/),(/1,1,1/) )
    
    L%x=L%x*PI/180.0_DP; L%y=L%y*PI/180.0_DP
    
    if (any(L%h(:,:,1:kmax-1)<0.0)) then
       call error_except('thickness of layer should not be less than 0')
    end if
    if (any(L%h(:,:,kmax)<SEIS_ZERO)) then
       print *, minloc(L%h(:,:,kmax))
       call error_except('thickness of lowest layer should larger than 0')
    end if
    
!-- density
    allocate(L%Dp(2,imax,jmax,kmax))
    allocate(L%Dp_poly_d(kmax))
#ifdef WITHQS
    allocate(L%Qs(2,imax,jmax,kmax))
    allocate(L%Qs_poly_d(kmax))
#endif
    call nfseis_varget(L%fnm,'rho',L%Dp,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
    call nfseis_varget(L%fnm,'rho_poly_d',L%Dp_poly_d,(/1/),(/kmax/),(/1/) )
#ifdef WITHQS
    call nfseis_varget(L%fnm,'Qs',L%Qs,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
    call nfseis_attget(L%fnm,'QsF0',L%QsF0)
    call nfseis_attget(L%fnm,'QsINF',L%QsINF)
    QsF0=L%QsF0; QsINF=L%QsINF
#endif
    
!-- stiffness from velocity parameters
    if (flag_stif_type==SIG_STIF_VEL) then
       allocate(L%Vp(2,imax,jmax,kmax))
       allocate(L%Vs(2,imax,jmax,kmax))
       allocate(L%Vp_poly_d(kmax))
       allocate(L%Vs_poly_d(kmax))
       if (flag_anis_type==SIG_ANIS_ISO) then
          call nfseis_varget(L%fnm,'Vp',L%Vp,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'Vs',L%Vs,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'Vp_poly_d',L%Vp_poly_d,(/1/),(/kmax/),(/1/) )
          call nfseis_varget(L%fnm,'Vs_poly_d',L%Vs_poly_d,(/1/),(/kmax/),(/1/) )
       end if
       if (flag_anis_type>=SIG_ANIS_VTI) then
          allocate(L%Vph(2,imax,jmax,kmax))
          allocate(L%Vsh(2,imax,jmax,kmax))
          allocate(L%eta(2,imax,jmax,kmax))
          allocate(L%Vph_poly_d(kmax))
          allocate(L%Vsh_poly_d(kmax))
          allocate(L%eta_poly_d(kmax))
          call nfseis_varget(L%fnm,'Vpv',L%Vp,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'Vsv',L%Vs,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'Vph',L%Vph,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'Vsh',L%Vsh,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'eta',L%eta,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'Vpv_poly_d',L%Vp_poly_d,(/1/),(/kmax/),(/1/) )
          call nfseis_varget(L%fnm,'Vsv_poly_d',L%Vs_poly_d,(/1/),(/kmax/),(/1/) )
          call nfseis_varget(L%fnm,'Vph_poly_d',L%Vph_poly_d,(/1/),(/kmax/),(/1/) )
          call nfseis_varget(L%fnm,'Vsh_poly_d',L%Vsh_poly_d,(/1/),(/kmax/),(/1/) )
          call nfseis_varget(L%fnm,'eta_poly_d',L%eta_poly_d,(/1/),(/kmax/),(/1/) )
       end if
       if (flag_anis_type>=SIG_ANIS_TTI) then
          allocate(L%theta(2,imax,jmax,kmax))
          allocate(L%phi(2,imax,jmax,kmax))
          allocate(L%theta_poly_d(kmax))
          allocate(L%phi_poly_d(kmax))
          call nfseis_varget(L%fnm,'TTI_theta',L%theta,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'TTI_phi',L%phi,(/1,1,1,1/),(/2,imax,jmax,kmax/),(/1,1,1,1/) )
          call nfseis_varget(L%fnm,'TTI_theta_poly_d',L%theta_poly_d,(/1/),(/kmax/),(/1/) )
          call nfseis_varget(L%fnm,'TTI_phi_poly_d',L%phi_poly_d,(/1/),(/kmax/),(/1/) )
          L%theta=L%theta*PI/180.0_DP
          L%phi=L%phi*PI/180.0_DP
       end if
       !-- TRICLINIC
       if (flag_anis_type==SIG_ANIS_TRICLINIC) then
          call error_except("no velocity equivalent for TRICLINIC medium")
       end if
    end if
    
!-- stiffness tensor
    if (flag_stif_type==SIG_STIF_TEN) then
       allocate(L%C(21,2,imax,jmax,kmax))
       allocate(L%C_poly_d(21,kmax))
       call nfseis_varget(L%fnm,'C',L%C,(/1,1,1,1,1/),(/21,2,imax,jmax,kmax/),(/1,1,1,1,1/) )
       call nfseis_varget(L%fnm,'C_poly_d',L%C_poly_d,(/1,1/),(/21,kmax/),(/1,1/) )
    end if
!-------------------------------------------------------------------------------
end subroutine composite_read
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine composite_destroy
!-------------------------------------------------------------------------------
    if (associated(LC%x  )         ) deallocate(LC%x  )
    if (associated(LC%y  )         ) deallocate(LC%y  )
    if (associated(LC%h  )         ) deallocate(LC%h  )
    if (associated(LC%Dp )         ) deallocate(LC%Dp )
    if (associated(LC%Dp_poly_d)   ) deallocate(LC%Dp_poly_d)
    if (associated(LC%Qs )         ) deallocate(LC%Qs )
    if (associated(LC%Qs_poly_d)   ) deallocate(LC%Qs_poly_d)
    if (associated(LC%Vp )         ) deallocate(LC%Vp )
    if (associated(LC%Vp_poly_d)   ) deallocate(LC%Vp_poly_d)
    if (associated(LC%Vs )         ) deallocate(LC%Vs )
    if (associated(LC%Vs_poly_d)   ) deallocate(LC%Vs_poly_d)
    if (associated(LC%Vph )        ) deallocate(LC%Vph )
    if (associated(LC%Vph_poly_d)  ) deallocate(LC%Vph_poly_d)
    if (associated(LC%Vsh )        ) deallocate(LC%Vsh )
    if (associated(LC%Vsh_poly_d)  ) deallocate(LC%Vsh_poly_d)
    if (associated(LC%eta )        ) deallocate(LC%eta )
    if (associated(LC%eta_poly_d)  ) deallocate(LC%eta_poly_d)
    if (associated(LC%theta )      ) deallocate(LC%theta )
    if (associated(LC%theta_poly_d)) deallocate(LC%theta_poly_d)
    if (associated(LC%phi )        ) deallocate(LC%phi )
    if (associated(LC%phi_poly_d)  ) deallocate(LC%phi_poly_d)
    if (associated(LC%C )          ) deallocate(LC%C )
    if (associated(LC%C_poly_d)    ) deallocate(LC%C_poly_d)
!-------------------------------------------------------------------------------
end subroutine composite_destroy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine composite_retrieve(x0,y0,d0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: x0,y0,d0,dz2
    real(SP),dimension(:),intent(out) :: Vp,Vs,Dp,Vph,Vsh,eta,theta,phi
    real(SP),dimension(:,:),intent(out) :: CMat
    
    real(SP) :: x1,y1
    real(SP) :: d,h,L1,L2
    
    integer :: n,m
    integer :: i1,i2,j1,j2,k1,k2,k0,n1,n2
    
    d=0.0
    call indx_locate_1d(x0,LC%x,i1,i2,x1)
    call indx_locate_1d(y0,LC%y,j1,j2,y1)
    do n=1,LC%nk
       h=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%h(i1:i2,j1:j2,n),2,2,x1,y1)
    
       if (h<=SEIS_ZERO) cycle
    
       d=d+h
       if (abs(d-d0)<=dz2/5.0) then
          if (n==LC%nk) call error_except('please increase the thickness of the lowest layer')
          n1=n
          do m=n+1,LC%nk
             h=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%h(i1:i2,j1:j2,m),2,2,x1,y1)
             if (h>SEIS_ZERO) then
                 n2=m
                 exit
             end if
          end do
          ! get value
          Dp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(2,i1:i2,j1:j2,n1),2,2,x1,y1)
          Dp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(1,i1:i2,j1:j2,n2),2,2,x1,y1)
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(2,i1:i2,j1:j2,n1),2,2,x1,y1)
             Vs(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(2,i1:i2,j1:j2,n1),2,2,x1,y1)
             Vp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(1,i1:i2,j1:j2,n2),2,2,x1,y1)
             Vs(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(1,i1:i2,j1:j2,n2),2,2,x1,y1)
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vph(2,i1:i2,j1:j2,n1),2,2,x1,y1)
                Vsh(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vsh(2,i1:i2,j1:j2,n1),2,2,x1,y1)
                eta(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%eta(2,i1:i2,j1:j2,n1),2,2,x1,y1)
                Vph(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vph(1,i1:i2,j1:j2,n2),2,2,x1,y1)
                Vsh(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vsh(1,i1:i2,j1:j2,n2),2,2,x1,y1)
                eta(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%eta(1,i1:i2,j1:j2,n2),2,2,x1,y1)
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%theta(2,i1:i2,j1:j2,n1),2,2,x1,y1)
                theta(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%theta(1,i1:i2,j1:j2,n2),2,2,x1,y1)
                phi(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%phi(2,i1:i2,j1:j2,n1),2,2,x1,y1)
                phi(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%phi(1,i1:i2,j1:j2,n2),2,2,x1,y1)
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(1,m)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%C(m,2,i1:i2,j1:j2,n1),2,2,x1,y1)
                CMat(2,m)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%C(m,1,i1:i2,j1:j2,n2),2,2,x1,y1)
             end do
          end if
          exit
       elseif (d>d0) then
          Dp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(1,i1:i2,j1:j2,n),2,2,x1,y1)
          Dp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Dp(2,i1:i2,j1:j2,n),2,2,x1,y1)
          Dp(:)=Dp(1)+(Dp(2)-Dp(1))/h**LC%Dp_poly_d(n)*(d0-(d-h))**LC%Dp_poly_d(n)
    
          if (flag_stif_type==SIG_STIF_VEL) then
             Vp(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(1,i1:i2,j1:j2,n),2,2,x1,y1)
             Vs(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(1,i1:i2,j1:j2,n),2,2,x1,y1)
             Vp(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vp(2,i1:i2,j1:j2,n),2,2,x1,y1)
             Vs(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vs(2,i1:i2,j1:j2,n),2,2,x1,y1)
             Vp(:)=Vp(1)+(Vp(2)-Vp(1))/h**LC%Vp_poly_d(n)*(d0-(d-h))**LC%Vp_poly_d(n)
             Vs(:)=Vs(1)+(Vs(2)-Vs(1))/h**LC%Vs_poly_d(n)*(d0-(d-h))**LC%Vs_poly_d(n)
             if (flag_anis_type>=SIG_ANIS_VTI) then
                Vph(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vph(1,i1:i2,j1:j2,n),2,2,x1,y1)
                Vph(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vph(2,i1:i2,j1:j2,n),2,2,x1,y1)
                Vph(:)=Vph(1)+(Vph(2)-Vph(1))/h**LC%Vph_poly_d(n)*(d0-(d-h))**LC%Vph_poly_d(n)
                Vsh(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vsh(1,i1:i2,j1:j2,n),2,2,x1,y1)
                Vsh(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%Vsh(2,i1:i2,j1:j2,n),2,2,x1,y1)
                Vsh(:)=Vsh(1)+(Vsh(2)-Vsh(1))/h**LC%Vsh_poly_d(n)*(d0-(d-h))**LC%Vsh_poly_d(n)
                eta(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%eta(1,i1:i2,j1:j2,n),2,2,x1,y1)
                eta(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%eta(2,i1:i2,j1:j2,n),2,2,x1,y1)
                eta(:)=eta(1)+(eta(2)-eta(1))/h**LC%eta_poly_d(n)*(d0-(d-h))**LC%eta_poly_d(n)
             end if
             if (flag_anis_type>=SIG_ANIS_TTI) then
                theta(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%theta(1,i1:i2,j1:j2,n),2,2,x1,y1)
                theta(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%theta(2,i1:i2,j1:j2,n),2,2,x1,y1)
                theta(:)=theta(1)+(theta(2)-theta(1))/h**LC%theta_poly_d(n)*(d0-(d-h))**LC%theta_poly_d(n)
                phi(1)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%phi(1,i1:i2,j1:j2,n),2,2,x1,y1)
                phi(2)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%phi(2,i1:i2,j1:j2,n),2,2,x1,y1)
                phi(:)=phi(1)+(phi(2)-phi(1))/h**LC%phi_poly_d(n)*(d0-(d-h))**LC%phi_poly_d(n)
             end if
          end if
          if (flag_stif_type==SIG_STIF_TEN) then
             do m=1,21
                CMat(1,m)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%C(m,1,i1:i2,j1:j2,n),2,2,x1,y1)
                CMat(2,m)=interp_2d(LC%x(i1:i2),LC%y(j1:j2),LC%C(m,2,i1:i2,j1:j2,n),2,2,x1,y1)
                CMat(:,m)=CMat(1,m)+(CMat(2,m)-CMat(1,m))/h**LC%C_poly_d(m,n)*(d0-(d-h))**LC%C_poly_d(m,n)
             end do
          end if
          exit
       elseif (n==LC%nk) then
          call error_except('please increase the thickness of the lowest layer')
       end if
    end do
!-------------------------------------------------------------------------------
end subroutine composite_retrieve
!-------------------------------------------------------------------------------

!===============================================================================
! -----------------------------  volume structure ----------------------------
!===============================================================================

!-------------------------------------------------------------------------------
subroutine volume_init(P)
!-------------------------------------------------------------------------------
    type(STRUCT_VOLUME) :: P
    character (len=SEIS_STRLEN) :: ctype
    integer :: imax,jmax,kmax
    
    if (trim(P%filetype)/='nc') then
       call error_except('only nc input accepted for volume type')
    end if
    
      call nfseis_diminfo(P%fnm,'theta',imax)
      call nfseis_diminfo(P%fnm,'phi',jmax)
      call nfseis_diminfo(P%fnm,'depth',kmax)
      allocate(P%gx(imax)); P%gx=0.0
      allocate(P%gy(jmax)); P%gy=0.0
      allocate(P%grz(kmax)); P%grz=0.0
      call nfseis_varget(P%fnm,'theta',P%gx,(/1/),(/imax/),(/1/))
      call nfseis_varget(P%fnm,'phi',P%gy,(/1/),(/jmax/),(/1/))
      call nfseis_varget(P%fnm,'depth2sealevel',P%grz,(/1/),(/kmax/),(/1/))
      call nfseis_attget(P%fnm,'sealevel',P%rz0)
    
      call nfseis_attget(P%fnm,'anisotrop_type',ctype)
      flag_anis_type=anis_name2id(trim(adjustl(ctype)))
      call nfseis_attget(P%fnm,'stiffness_type',ctype)
      flag_stif_type=stif_name2id(trim(adjustl(ctype)))
    
      !== convert degree to radian ==
      P%gx=P%gx*PI/180.0_DP; P%gy=P%gy*PI/180.0_DP
      !== save dimension information ==
      P%imax=imax; P%jmax=jmax; P%kmax=kmax
    
      !== pre-allocate P ==
      !!== suppose the dencest input is sampled as the effective vol-integrate
      !!== point, which is equal to nx*NGRIDX, etc. Allocate minimal required
      !!== grid for reading into this thread
      !imax=min(nx*NGRIDX,imax)
      !jmax=min(ny*NGRIDY,jmax)
      !kmax=min(nz*NGRIDZ,kmax)
      !call volume_alloc(P,imax,jmax,kmax)
!-------------------------------------------------------------------------------
end subroutine volume_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine volume_alloc(P,imax,jmax,kmax)
!-------------------------------------------------------------------------------
    type(STRUCT_VOLUME) :: P
    integer,intent(in) :: imax,jmax,kmax
    
    if (flag_stif_type==SIG_STIF_VEL) then
       if (associated(P%Vp)) then
       if (size(P%Vp,1)<imax .or. size(P%Vp,2)<jmax .or. size(P%Vp,3)<kmax) then
          deallocate(P%x);deallocate(P%y);deallocate(P%rz);
          deallocate(P%Vp); deallocate(P%Vs); deallocate(P%Dp)
          if (flag_anis_type>=SIG_ANIS_VTI) then
             deallocate(P%Vph); deallocate(P%Vsh); deallocate(P%eta)
          end if
          if (flag_anis_type>=SIG_ANIS_TTI) then
             deallocate(P%theta); deallocate(P%phi)
          end if
       end if
       end if
       
       if (.not. associated(P%Vp)) then
          allocate(P%x(imax)); P%x=0.0
          allocate(P%y(jmax)); P%y=0.0
          allocate(P%rz(kmax)); P%rz=0.0
          allocate(P%Vp(imax,jmax,kmax)); P%Vp=0.0
          allocate(P%Vs(imax,jmax,kmax)); P%Vs=0.0
          allocate(P%Dp(imax,jmax,kmax)); P%Dp=0.0
          if (flag_anis_type>=SIG_ANIS_VTI) then
             allocate(P%Vph(imax,jmax,kmax)); P%Vph=0.0
             allocate(P%Vsh(imax,jmax,kmax)); P%Vsh=0.0
             allocate(P%eta(imax,jmax,kmax)); P%eta=0.0
          end if
          if (flag_anis_type>=SIG_ANIS_TTI) then
             allocate(P%theta(imax,jmax,kmax)); P%theta=0.0
             allocate(P%phi(imax,jmax,kmax)); P%phi=0.0
          end if
       end if
    end if
    
    if (flag_stif_type==SIG_STIF_TEN) then
       if (associated(P%C)) then
       if (size(P%C,2)<imax .or. size(P%C,3)<jmax .or. size(P%C,4)<kmax) then
          deallocate(P%x);deallocate(P%y);deallocate(P%rz);
          deallocate(P%Dp)
          deallocate(P%C)
       end if
       end if
       
       if (.not. associated(P%Vp)) then
          allocate(P%x(imax)); P%x=0.0
          allocate(P%y(jmax)); P%y=0.0
          allocate(P%rz(kmax)); P%rz=0.0
          allocate(P%Dp(imax,jmax,kmax)); P%Dp=0.0
          allocate(P%C(21,imax,jmax,kmax)); P%C=0.0
       end if
    end if
    
    P%ni=imax;P%nj=jmax;P%nk=kmax
!-------------------------------------------------------------------------------
end subroutine volume_alloc
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine volume_read(P)
!-------------------------------------------------------------------------------
    type(STRUCT_VOLUME) :: P
    integer :: i1,i2,j1,j2,k1,k2,imax,jmax,kmax
    integer :: indx(1)
    indx=maxloc(P%gx,P%gx<=max(x(nx1),P%gx(1)));i1=max(indx(1),1)
    indx=minloc(P%gx,P%gx>=min(x(nx2),P%gx(P%imax)));i2=indx(1)
    indx=maxloc(P%gy,P%gy<=max(y(ny1),P%gy(1)));j1=max(indx(1),1)
    indx=minloc(P%gy,P%gy>=min(y(ny2),P%gy(P%jmax)));j2=indx(1)
    k1=1; k2=P%kmax
    imax=i2-i1+1; jmax=j2-j1+1; kmax=k2-k1+1
    call volume_alloc(P,imax,jmax,kmax)
    
    P%x(1:imax)=P%gx(i1:i2); P%y(1:jmax)=P%gy(j1:j2)
    P%rz(1:kmax)=P%grz(k1:k2);
    
!-- rho
    call nfseis_varget( P%fnm,'rho',P%Dp(1:imax,1:jmax,1:kmax), &
         (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
    
!-- stiffness from velocity parameters
    if (flag_stif_type==SIG_STIF_VEL) then
       if (flag_anis_type==SIG_ANIS_ISO) then
          call nfseis_varget( P%fnm,'Vp',P%Vp(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          call nfseis_varget( P%fnm,'Vs',P%Vs(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
       end if
       if (flag_anis_type>=SIG_ANIS_VTI) then
          call nfseis_varget( P%fnm,'Vpv',P%Vp(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          call nfseis_varget( P%fnm,'Vsv',P%Vs(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          call nfseis_varget( P%fnm,'Vph',P%Vph(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          call nfseis_varget( P%fnm,'Vsh',P%Vsh(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          call nfseis_varget( P%fnm,'eta',P%eta(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
       end if
       if (flag_anis_type>=SIG_ANIS_TTI) then
          call nfseis_varget( P%fnm,'TTI_theta',P%theta(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          call nfseis_varget( P%fnm,'TTI_phi',P%phi(1:imax,1:jmax,1:kmax), &
               (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
          P%theta=P%theta*PI/180.0_DP
          P%phi=P%phi*PI/180.0_DP
       end if
    end if
!-- stiffness tensor
    if (flag_stif_type==SIG_STIF_TEN) then
       call nfseis_varget( P%fnm,'C',P%C(1:21,1:imax,1:jmax,1:kmax), &
            (/1,i1,j1,k1/),(/21,imax,jmax,kmax/),(/1,1,1,1/) )
    end if
!-------------------------------------------------------------------------------
end subroutine volume_read
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine volume_destroy
!-------------------------------------------------------------------------------
    if (associated(BV%gx  )  ) deallocate(BV%gx  )
    if (associated(BV%gy  )  ) deallocate(BV%gy  )
    if (associated(BV%grz )  ) deallocate(BV%grz )
    if (associated(BV%x  )  ) deallocate(BV%x  )
    if (associated(BV%y  )  ) deallocate(BV%y  )
    if (associated(BV%rz )  ) deallocate(BV%rz )
    if (associated(BV%Dp )  ) deallocate(BV%Dp )
    if (associated(BV%Vp )  ) deallocate(BV%Vp )
    if (associated(BV%Vs )  ) deallocate(BV%Vs )
    if (associated(BV%Vph)  ) deallocate(BV%Vph)
    if (associated(BV%Vsh)  ) deallocate(BV%Vsh)
    if (associated(BV%eta)  ) deallocate(BV%eta)
    if (associated(BV%theta))   deallocate(BV%theta)
    if (associated(BV%phi)  ) deallocate(BV%phi)
    if (associated(BV%C)    ) deallocate(BV%C)
!-------------------------------------------------------------------------------
end subroutine volume_destroy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine volume_retrieve(x0,y0,z0,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: x0,y0,z0
    real(SP),dimension(:),intent(out) :: Vp,Vs,Dp,Vph,Vsh,eta,theta,phi
    real(SP),dimension(:,:),intent(out) :: CMat
    
    real(SP) :: x1,y1,z1
    
    integer :: i1,i2,j1,j2,k1,k2
    integer :: m
    
    call indx_locate_1d(x0,BV%x(1:BV%ni),i1,i2,x1)
    call indx_locate_1d(y0,BV%y(1:BV%nj),j1,j2,y1)
    call indx_locate_1d(BV%rz0-z0,BV%rz(1:BV%nk),k1,k2,z1)
    
    Dp(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
          BV%Dp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
    if (flag_stif_type==SIG_STIF_VEL) then
       Vp(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
             BV%Vp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
       Vs(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
             BV%Vs(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
       if (flag_anis_type>=SIG_ANIS_VTI) then
          Vph(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
                 BV%Vph(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
          Vsh(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
                 BV%Vsh(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
          eta(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
                 BV%eta(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
       end if
       if (flag_anis_type>=SIG_ANIS_TTI) then
          theta(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
                 BV%theta(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
          phi(:)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
                 BV%phi(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
       end if
    end if
    if (flag_stif_type==SIG_STIF_TEN) then
       do m=1,21
          CMat(:,m)=interp_3d(BV%x(i1:i2),BV%y(j1:j2),BV%rz(k1:k2), &
                    BV%C(m,i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
       end do
    end if
!-------------------------------------------------------------------------------
end subroutine volume_retrieve
!-------------------------------------------------------------------------------

!===============================================================================
! ---------------------------  verpoly structure ----------------------------
!===============================================================================

!-------------------------------------------------------------------------------
subroutine verpoly_init(P)
!-------------------------------------------------------------------------------
    type(STRUCT_VERPOLY) :: P
    integer :: imax,jmax,kmax,nmax
    
    if (trim(P%filetype)/='nc') then
       call error_except('only nc input accepted for verpoly type')
    end if
    
    call nfseis_diminfo(P%fnm,'theta',imax)
    call nfseis_diminfo(P%fnm,'phi',jmax)
    call nfseis_diminfo(P%fnm,'depth',kmax)
    call nfseis_diminfo(P%fnm,'polynomial',nmax)
    allocate(P%gx(imax)); P%gx=0.0
    allocate(P%gy(jmax)); P%gy=0.0
    allocate(P%Vp_poly_d(nmax)); P%Vp_poly_d=0.0
    allocate(P%Vs_poly_d(nmax)); P%Vs_poly_d=0.0
    allocate(P%Dp_poly_d(nmax)); P%Dp_poly_d=0.0
    call nfseis_varget(P%fnm,'theta',P%gx,(/1/),(/imax/),(/1/))
    call nfseis_varget(P%fnm,'phi',P%gy,(/1/),(/jmax/),(/1/))
    call nfseis_varget(P%fnm,'Vp_poly_d',P%Vp_poly_d,(/1/),(/nmax/),(/1/))
    call nfseis_varget(P%fnm,'Vs_poly_d',P%Vs_poly_d,(/1/),(/nmax/),(/1/))
    call nfseis_varget(P%fnm,'rho_poly_d',P%Dp_poly_d,(/1/),(/nmax/),(/1/))
    call nfseis_attget(P%fnm,'sealevel',P%rz0)
    P%gx=P%gx*PI/180.0_DP; P%gy=P%gy*PI/180.0_DP
    P%imax=imax; P%jmax=jmax; P%kmax=kmax; P%nmax=nmax
    imax=min(nx*max(2*NGRIDX,1),imax)
    jmax=min(ny*max(2*NGRIDY,1),jmax)
    kmax=min(nz*max(2*NGRIDZ,1),kmax)
    call verpoly_alloc(P,imax,jmax,kmax)
!-------------------------------------------------------------------------------
end subroutine verpoly_init
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine verpoly_alloc(P,imax,jmax,kmax)
!-------------------------------------------------------------------------------
    type(STRUCT_VERPOLY) :: P
    integer,intent(in) :: imax,jmax,kmax
    
    if (associated(P%rz)) then
    if (size(P%rz,1)<imax .or. size(P%rz,2)<jmax .or. size(P%rz,3)<kmax) then
       deallocate(P%x);deallocate(P%y);deallocate(P%rz);
       deallocate(P%Vp_poly_c); deallocate(P%Vs_poly_c); deallocate(P%Dp_poly_c)
    end if
    end if
    
    P%ni=imax;P%nj=jmax;P%nk=kmax
    if (.not. associated(P%rz)) then
       allocate(P%x(imax)); P%x=0.0
       allocate(P%y(jmax)); P%y=0.0
       allocate(P%rz(imax,jmax,kmax)); P%rz=0.0
       allocate(P%Vp_poly_c(imax,jmax,kmax,P%nmax)); P%Vp_poly_c=0.0
       allocate(P%Vs_poly_c(imax,jmax,kmax,P%nmax)); P%Vs_poly_c=0.0
       allocate(P%Dp_poly_c(imax,jmax,kmax,P%nmax)); P%Dp_poly_c=0.0
    end if
!-------------------------------------------------------------------------------
end subroutine verpoly_alloc
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine verpoly_read(P)
!-------------------------------------------------------------------------------
    type(STRUCT_VERPOLY) :: P
    integer :: i1,i2,j1,j2,k1,k2,imax,jmax,kmax,nmax
    integer :: indx(1)
    indx=maxloc(P%gx,P%gx<=max(x(nx1),P%gx(1)));i1=max(indx(1),1)
    indx=minloc(P%gx,P%gx>=min(x(nx2),P%gx(P%imax)));i2=indx(1)
    indx=maxloc(P%gy,P%gy<=max(y(ny1),P%gy(1)));j1=max(indx(1),1)
    indx=minloc(P%gy,P%gy>=min(y(ny2),P%gy(P%jmax)));j2=indx(1)
    k1=1; k2=P%kmax
    imax=i2-i1+1; jmax=j2-j1+1; kmax=k2-k1+1; nmax=P%nmax
    call verpoly_alloc(P,imax,jmax,kmax)
    
    P%x(1:imax)=P%gx(i1:i2); P%y(1:jmax)=P%gy(j1:j2)
    call nfseis_varget( P%fnm,'depth2sealevel',P%rz(1:imax,1:jmax,1:kmax), &
         (/i1,j1,k1/),(/imax,jmax,kmax/),(/1,1,1/) )
    call nfseis_varget( P%fnm,'Vp_poly_c',P%Vp_poly_c(1:imax,1:jmax,1:kmax,:), &
         (/i1,j1,k1,1/),(/imax,jmax,kmax,nmax/),(/1,1,1,1/) )
    call nfseis_varget( P%fnm,'Vs_poly_c',P%Vs_poly_c(1:imax,1:jmax,1:kmax,:), &
         (/i1,j1,k1,1/),(/imax,jmax,kmax,nmax/),(/1,1,1,1/) )
    call nfseis_varget( P%fnm,'rho_poly_c',P%Dp_poly_c(1:imax,1:jmax,1:kmax,:), &
         (/i1,j1,k1,1/),(/imax,jmax,kmax,nmax/),(/1,1,1,1/) )
!-------------------------------------------------------------------------------
end subroutine verpoly_read
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine verpoly_destroy
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
end subroutine verpoly_destroy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine verpoly_retrieve(x0,y0,z0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: x0,y0,z0,dz2
    real(SP),dimension(:),intent(out) :: Vp,Vs,Dp,Vph,Vsh,eta,theta,phi
    real(SP),dimension(:,:),intent(out) :: CMat
    
    real(SP) :: x1,y1,z1,z2
    
    integer :: indx(1)
    integer :: k1,k2,wi,wj,wk
    integer :: m
    
    indx=minloc(abs(BP%x(1:BP%ni)-x0)); wi=indx(1)
    indx=minloc(abs(BP%y(1:BP%nj)-y0)); wj=indx(1)
    !if (abs(BP%x(wi)-x0)+abs(BP%y(wj)-y0)>SEIS_ZERO) then
    !   print *, 'error'
    !   stop 1
    !else
    indx=minloc(abs(BP%rz(wi,wj,1:BP%nk)-(BP%rz0-z0))); wk=indx(1)
    if (abs(BP%rz(wi,wj,wk)-(BP%rz0-z0))<=dz2/5.0  &
       .or. (BP%rz0-z0)<BP%rz(wi,wj,1)) then
       k1=max(1,wk-1); k2=wk
       z1=BP%rz(wi,wj,k2)-BP%rz(wi,wj,k1); z2=0
       Vp(1)=eval_poly(BP%Vp_poly_c(wi,wj,k1,:),BP%Vp_poly_d,z1)
       Vp(2)=eval_poly(BP%Vp_poly_c(wi,wj,k2,:),BP%Vp_poly_d,z2)
       Vs(1)=eval_poly(BP%Vs_poly_c(wi,wj,k1,:),BP%Vs_poly_d,z1)
       Vs(2)=eval_poly(BP%Vs_poly_c(wi,wj,k2,:),BP%Vs_poly_d,z2)
       Dp(1)=eval_poly(BP%Dp_poly_c(wi,wj,k1,:),BP%Dp_poly_d,z1)
       Dp(2)=eval_poly(BP%Dp_poly_c(wi,wj,k2,:),BP%Dp_poly_d,z2)
    else
       if (BP%rz0-z0<BP%rz(wi,wj,wk)) wk=wk-1
       z1=(BP%rz0-z0)-BP%rz(wi,wj,wk)
       Vp(:)=eval_poly(BP%Vp_poly_c(wi,wj,wk,:),BP%Vp_poly_d,z1)
       Vs(:)=eval_poly(BP%Vs_poly_c(wi,wj,wk,:),BP%Vs_poly_d,z1)
       Dp(:)=eval_poly(BP%Dp_poly_c(wi,wj,wk,:),BP%Dp_poly_d,z1)
    end if
!-------------------------------------------------------------------------------
end subroutine verpoly_retrieve
!-------------------------------------------------------------------------------

!===============================================================================
subroutine effmedia_eval
!===============================================================================
!-- This subroutine performs effective parameter volume integration   ==
!-- by sampling medium parameters at several points in the cell, then ==
!-- summation and average ==

    real(SP),dimension(NGRIDX) :: xvec
    real(SP),dimension(NGRIDY) :: yvec
    real(SP),dimension(NGRIDZ) :: zvec
    real(SP) :: dx1,dx2,dy1,dy2,dz1,dz2
    
    integer :: ncrit                            != criteria to average stiffness
    integer,dimension(21) :: nsamped            != summation number
    real(SP),dimension(21) :: vsamped           != summation value
    real(SP),dimension(6,6) :: TTIC,TTIM
    real(SP) :: rho0
    
    !real(SP),dimension(2) :: lam, miu
    real(SP),dimension(2,21) :: CMat
    real(SP),dimension(2) :: Vp,Vs,Dp,dVp,dVs,dDp
    real(SP),dimension(2) :: Vph,Vsh,eta
    real(SP),dimension(2) :: theta,phi
    
    real(SP) :: x0,y0,z0,d0
    !real(SP) :: x1,y1,z1,z2
    integer :: nsamp
    integer :: i,j,k,mi,mj,mk,n,m
    !integer :: i1,i2,j1,j2,k1,k2,k0,n1,n2
    !integer :: indx(1)
    
    !layered
    !real(SP) :: h,d
    !real(SP) :: L1,L2
    !verpoly
    !integer :: wi,wj,wk
    
    !== total sampling points for volume integration ==
    nsamp=NGRIDX*NGRIDY*NGRIDZ

!-------------------------------------------------------------------------------
!-- loop all grid points --
!-------------------------------------------------------------------------------
do k=nk1,nk2
#ifdef MediaMPI
   if (masternode)  &
#endif
       print *, ' k=',k-nk1+1, ' of',nk
do j=nj1,nj2
do i=ni1,ni2

   !== init count and value of summation to zero ==
   nsamped=0; vsamped=0.0_SP; rho0=0.0;

   !==== the sampling position for one cell ====
   !-- half length of each sampled interval --
   dx1=(x(i)-x(i-1))/2.0/NGRIDX; dx2=(x(i+1)-x(i))/2.0/NGRIDX
   !-- x < x(i) --
   do n=1,(NGRIDX+1)/2
      xvec(n)=(x(i)+x(i-1))/2.0+dx1*(2*n-1)
   end do
   !-- x > x(i)
   do n=(NGRIDX+1)/2+1,NGRIDX
      xvec(n)=(x(i)+x(i+1))/2.0-dx2*(2*(NGRIDX-n+1)-1)
   end do
   !-- simliar for y --
   dy1=(y(j)-y(j-1))/2.0/NGRIDY; dy2=(y(j+1)-y(j))/2.0/NGRIDY
   do n=1,(NGRIDY+1)/2
      yvec(n)=(y(j)+y(j-1))/2.0+dy1*(2*n-1)
   end do
   do n=(NGRIDY+1)/2+1,NGRIDY
      yvec(n)=(y(j)+y(j+1))/2.0-dy2*(2*(NGRIDY-n+1)-1)
   end do
   !-- simliar for z --
   dz1=(z(k)-z(k-1))/2.0/NGRIDZ; dz2=(z(k+1)-z(k))/2.0/NGRIDZ
   do n=1,(NGRIDZ+1)/2
      zvec(n)=(z(k)+z(k-1))/2.0+dz1*(2*n-1)
   end do
   do n=(NGRIDZ+1)/2+1,NGRIDZ
      zvec(n)=(z(k)+z(k+1))/2.0-dz2*(2*(NGRIDZ-n+1)-1)
   end do

!-------------------------------------------------------------------------------
!-- loop all sampling position in one cell --
!-------------------------------------------------------------------------------
   do mk=1,NGRIDZ
   do mj=1,NGRIDY
   do mi=1,NGRIDX
   
      CMat=0.0_SP
      x0=xvec(mi); y0=yvec(mj); z0=zvec(mk);
      d0=max(ztopo-z0,0.0);
   
!--  background model --
   !-- cart1d
   if (L1D%yes) then
      call layer1d_retrieve(d0,Vp,Vs,Dp,Vph,Vsh,eta)
   end if
   !-- interface
   if (LA%yes) then
      call interface_retrieve(x0,y0,z0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
   end if
   !-- layered
   if (L3D%yes) then
      call layered_retrieve(x0,y0,d0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
   end if
   !-- composite
   if (LC%yes) then
      call composite_retrieve(x0,y0,d0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
   end if
   !-- volume
   if (BV%yes) then
      call volume_retrieve(x0,y0,z0,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
   end if
   !-- verpoly
   if (BP%yes) then
      call verpoly_retrieve(x0,y0,z0,dz2,Vp,Vs,Dp,Vph,Vsh,eta,theta,phi,CMat)
   end if
   
!-- convert to stiffness tensor if model given in velocity-like --
   if (flag_stif_type==SIG_STIF_VEL) then
      if (flag_anis_type==SIG_ANIS_ISO) then
         CMat(:,21)=Dp*Vs**2  !mu
         CMat(:, 3)=Dp*Vp**2-2.0*CMat(:,21)  !lambda
#if defined AnisVTI || defined AnisGene
         !== for input ISO media in compiled code with anisotropy ==
         CMat(:, 1)=Dp*Vp**2  !A
         CMat(:,12)=Dp*Vp**2  !C
         CMat(:,16)=Dp*Vs**2  !L
#endif
      end if
      if (flag_anis_type>=SIG_ANIS_VTI) then
         CMat(:, 1)=Dp*Vph**2  !A
         CMat(:,12)=Dp*Vp **2  !C
         CMat(:,16)=Dp*Vs **2  !L
         CMat(:,21)=Dp*Vsh**2  !N
         CMat(:, 3)=eta*(CMat(:,1)-2.0_SP*CMat(:,16)) !F
      end if
   end if
   
   !-- rotate tensor for TTI media --
   if (flag_anis_type==SIG_ANIS_TTI) then
      do m=1,2
         TTIC=0.0_SP
         TTIC(1,1)=CMat(m, 1)
         TTIC(1,2)=CMat(m, 1)-2.0_SP*CMat(m, 21)
         TTIC(1,3)=CMat(m, 3)
         TTIC(2,2)=CMat(m, 1)
         TTIC(2,3)=CMat(m, 3)
         TTIC(3,3)=CMat(m,12)
         TTIC(4,4)=CMat(m,16)
         TTIC(5,5)=CMat(m,16)
         TTIC(6,6)=CMat(m,21)
         TTIC(2,1)=TTIC(1,2)
         TTIC(3,1)=TTIC(1,3)
         TTIC(3,2)=TTIC(2,3)
         TTIM=build_Bond(theta(m),phi(m))
         TTIC=matmul(matmul(TTIM,TTIC),transpose(TTIM))
         CMat(m, 1)=TTIC(1,1)
         CMat(m, 2)=TTIC(1,2)
         CMat(m, 3)=TTIC(1,3)
         CMat(m, 4)=TTIC(1,4)
         CMat(m, 5)=TTIC(1,5)
         CMat(m, 6)=TTIC(1,6)
         CMat(m, 7)=TTIC(2,2)
         CMat(m, 8)=TTIC(2,3)
         CMat(m, 9)=TTIC(2,4)
         CMat(m,10)=TTIC(2,5)
         CMat(m,11)=TTIC(2,6)
         CMat(m,12)=TTIC(3,3)
         CMat(m,13)=TTIC(3,4)
         CMat(m,14)=TTIC(3,5)
         CMat(m,15)=TTIC(3,6)
         CMat(m,16)=TTIC(4,4)
         CMat(m,17)=TTIC(4,5)
         CMat(m,18)=TTIC(4,6)
         CMat(m,19)=TTIC(5,5)
         CMat(m,20)=TTIC(5,6)
         CMat(m,21)=TTIC(6,6)
      end do
   end if
   
!------------- accumulate -------------
      rho0=rho0+0.5*(Dp(1)+Dp(2))
      do n=1,21
         !== only harmonically summate non-negative stiffness tensor? ==
         if (CMat(1,n)>=SEIS_ZERO .and. CMat(2,n)>=SEIS_ZERO) then
            vsamped(n)=vsamped(n)+0.5_SP/CMat(1,n)+0.5_SP/CMat(2,n)
            nsamped(n)=nsamped(n)+1
         end if
      end do
   
   end do !mi
   end do !mj
   end do !mk
!-------------------------------------------------------------------------------
   
   ! 1  2  3  4  5  6  #1
   !    7  8  9 10 11  #2
   !      12 13 14 15  #3
   !         16 17 18  #4
   !            19 20  #5
   !               21  #6
   
    !== arithmetic average of density ==
    rho(i,j,k)=rho0/nsamp

    !== harmonic average of stiffness tensor if non-zero value ==
    !== larger then a criteria value of ncrit ==
    ncrit=nsamp/2
    !ncrit=0
                           C13(i,j,k)=0.0_SP
    if (nsamped( 3)>ncrit) C13(i,j,k)=nsamped( 3)/vsamped( 3)
                           C66(i,j,k)=0.0_SP
    if (nsamped(21)>ncrit) C66(i,j,k)=nsamped(21)/vsamped(21)
#if defined AnisVTI || defined AnisGene
                           C11(i,j,k)=0.0_SP
    if (nsamped( 1)>ncrit) C11(i,j,k)=nsamped( 1)/vsamped( 1)
                           C33(i,j,k)=0.0_SP
    if (nsamped(12)>ncrit) C33(i,j,k)=nsamped(12)/vsamped(12)
                           C44(i,j,k)=0.0_SP
    if (nsamped(16)>ncrit) C44(i,j,k)=nsamped(16)/vsamped(16)
#endif
#ifdef AnisGene
                           C12(i,j,k)=0.0_SP
    if (nsamped( 2)>ncrit) C12(i,j,k)=nsamped( 2)/vsamped( 2)
                           C14(i,j,k)=0.0_SP
    if (nsamped( 4)>ncrit) C14(i,j,k)=nsamped( 4)/vsamped( 4)
                           C15(i,j,k)=0.0_SP
    if (nsamped( 5)>ncrit) C15(i,j,k)=nsamped( 5)/vsamped( 5)
                           C16(i,j,k)=0.0_SP
    if (nsamped( 6)>ncrit) C16(i,j,k)=nsamped( 6)/vsamped( 6)
                           C22(i,j,k)=0.0_SP
    if (nsamped( 7)>ncrit) C22(i,j,k)=nsamped( 7)/vsamped( 7)
                           C23(i,j,k)=0.0_SP
    if (nsamped( 8)>ncrit) C23(i,j,k)=nsamped( 8)/vsamped( 8)
                           C24(i,j,k)=0.0_SP
    if (nsamped( 9)>ncrit) C24(i,j,k)=nsamped( 9)/vsamped( 9)
                           C25(i,j,k)=0.0_SP
    if (nsamped(10)>ncrit) C25(i,j,k)=nsamped(10)/vsamped(10)
                           C26(i,j,k)=0.0_SP
    if (nsamped(11)>ncrit) C26(i,j,k)=nsamped(11)/vsamped(11)
                           C34(i,j,k)=0.0_SP
    if (nsamped(13)>ncrit) C34(i,j,k)=nsamped(13)/vsamped(13)
                           C35(i,j,k)=0.0_SP
    if (nsamped(14)>ncrit) C35(i,j,k)=nsamped(14)/vsamped(14)
                           C36(i,j,k)=0.0_SP
    if (nsamped(15)>ncrit) C36(i,j,k)=nsamped(15)/vsamped(15)
                           C45(i,j,k)=0.0_SP
    if (nsamped(17)>ncrit) C45(i,j,k)=nsamped(17)/vsamped(17)
                           C46(i,j,k)=0.0_SP
    if (nsamped(18)>ncrit) C46(i,j,k)=nsamped(18)/vsamped(18)
                           C55(i,j,k)=0.0_SP
    if (nsamped(19)>ncrit) C55(i,j,k)=nsamped(19)/vsamped(19)
                           C56(i,j,k)=0.0_SP
    if (nsamped(20)>ncrit) C56(i,j,k)=nsamped(20)/vsamped(20)
#endif

end do !i
end do !j
end do !k
!===============================================================================
end subroutine effmedia_eval
!===============================================================================

!! -------------   perturbed model ---------------------------
!! volume
!if (PV%yes) then
!   call indx_locate_1d(x0,PV%x(1:PV%ni),i1,i2,x1)
!   call indx_locate_1d(y0,PV%y(1:PV%nj),j1,j2,y1)
!   call indx_locate_1d(PV%rz0-z0,PV%rz(1:PV%nk),k1,k2,z1)
!   dVp=interp_3d(PV%x(i1:i2),PV%y(j1:j2),PV%rz(k1:k2), &
!         PV%Vp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
!   dVs=interp_3d(PV%x(i1:i2),PV%y(j1:j2),PV%rz(k1:k2), &
!         PV%Vs(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
!   dDp=interp_3d(PV%x(i1:i2),PV%y(j1:j2),PV%rz(k1:k2), &
!         PV%Dp(i1:i2,j1:j2,k1:k2),2,2,2,x1,y1,z1)
!   Vp=Vp*(1.0+dVp)
!   Vs=Vs*(1.0+dVs)
!   Dp=Dp*(1.0+dDp)
!end if
!! verpoly
!if (PP%yes) then
!   indx=minloc(abs(PP%x(1:PP%ni)-x0)); wi=indx(1)
!   indx=minloc(abs(PP%y(1:PP%nj)-y0)); wj=indx(1)
!   !if (abs(PP%x(wi)-x0)+abs(PP%y(wj)-y0)>SEIS_ZERO) then
!   !   print *, 'error'
!   !   stop 1
!   !else
!   indx=minloc(abs(PP%rz(wi,wj,1:PP%nk)-(PP%rz0-z0))); wk=indx(1)
!   if (abs(PP%rz(wi,wj,wk)-(PP%rz0-z0))<=dz2/5.0  &
!      .or. (PP%rz0-z0)<PP%rz(wi,wj,1)) then
!      k1=max(1,wk-1); k2=wk
!      z1=PP%rz(wi,wj,k2)-PP%rz(wi,wj,k1); z2=0
!      dVp(1)=eval_poly(PP%Vp_poly_c(wi,wj,k1,:),PP%Vp_poly_d,z1)
!      dVp(2)=eval_poly(PP%Vp_poly_c(wi,wj,k2,:),PP%Vp_poly_d,z2)
!      dVs(1)=eval_poly(PP%Vs_poly_c(wi,wj,k1,:),PP%Vs_poly_d,z1)
!      dVs(2)=eval_poly(PP%Vs_poly_c(wi,wj,k2,:),PP%Vs_poly_d,z2)
!      dDp(1)=eval_poly(PP%Dp_poly_c(wi,wj,k1,:),PP%Dp_poly_d,z1)
!      dDp(2)=eval_poly(PP%Dp_poly_c(wi,wj,k2,:),PP%Dp_poly_d,z2)
!   else
!      if (PP%rz0-z0<PP%rz(wi,wj,wk)) wk=wk-1
!      z1=(PP%rz0-z0)-PP%rz(wi,wj,wk)
!      dVp(:)=eval_poly(PP%Vp_poly_c(wi,wj,wk,:),PP%Vp_poly_d,z1)
!      dVs(:)=eval_poly(PP%Vs_poly_c(wi,wj,wk,:),PP%Vs_poly_d,z1)
!      dDp(:)=eval_poly(PP%Dp_poly_c(wi,wj,wk,:),PP%Dp_poly_d,z1)
!   end if
!   Vp=Vp*(1.0+dVp)
!   Vs=Vs*(1.0+dVs)
!   Dp=Dp*(1.0+dDp)
!end if

!-------------------------------------------------------------------------------
function build_Bond(theta,phi) result(M)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: theta,phi
    real(SP),dimension(6,6) :: M
    real(SP),dimension(3,3) :: A
    A(1,:)=(/ cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta) /)
    A(2,:)=(/ -cos(phi)          , cos(phi)           , 0.0_SP      /)
    A(3,:)=(/ cos(phi)*sin(theta), sin(phi)*sin(theta),  cos(theta) /)
    M(1,:)=(/ A(1,1)**2, A(1,2)**2, A(1,3)**2, 2.0*A(1,2)*A(1,3), 2.0*A(1,3)*A(1,1), 2*A(1,1)*A(1,2) /)
    M(2,:)=(/ A(2,1)**2, A(2,2)**2, A(2,3)**2, 2.0*A(2,2)*A(2,3), 2.0*A(2,3)*A(2,1), 2*A(2,1)*A(2,2) /)
    M(3,:)=(/ A(3,1)**2, A(3,2)**2, A(3,3)**2, 2.0*A(3,2)*A(3,3), 2.0*A(3,3)*A(3,1), 2*A(3,1)*A(3,2) /)
    M(4,:)=(/ A(2,1)*A(3,1),A(2,2)*A(3,2),A(2,3)*A(3,3),A(2,2)*A(3,3)+A(2,3)*A(3,2),A(2,1)*A(3,3)+A(2,3)*A(3,1),A(2,2)*A(3,1)+A(2,1)*A(3,2)/)
    M(5,:)=(/ A(3,1)*A(1,1),A(3,2)*A(1,2),A(3,3)*A(1,3),A(1,2)*A(3,3)+A(1,3)*A(3,2),A(1,3)*A(3,1)+A(1,1)*A(3,3),A(1,1)*A(3,2)+A(1,2)*A(3,1)/)
    M(6,:)=(/ A(1,1)*A(2,1),A(1,2)*A(2,2),A(1,3)*A(2,3),A(1,2)*A(2,3)+A(1,3)*A(2,2),A(1,3)*A(2,1)+A(1,1)*A(2,3),A(1,1)*A(2,2)+A(1,2)*A(2,1)/)
!-------------------------------------------------------------------------------
end function build_Bond
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function eval_poly(c,d,x) result(f)
!-------------------------------------------------------------------------------
    real(SP),dimension(:),intent(in) :: c,d
    real(SP),intent(in) :: x
    real(SP) :: f
    integer :: nmax,n
    nmax=size(d)
    !f=pp(1)*x0**3+pp(2)*x0**2+pp(3)*x0+pp(4)
    f=0.0
    do n=1,nmax
       f=f+c(n) * (x**d(n))
    end do
!-------------------------------------------------------------------------------
end function eval_poly
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine media_extend
!-------------------------------------------------------------------------------
    call extend_equal(rho)
    call extend_equal(C66)
    call extend_equal(C13)
!-------------------------------------------------------------------------------
end subroutine media_extend
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine extend_equal(w)
!-------------------------------------------------------------------------------
    real(SP),dimension(nx1:nx2,ny1:ny2,nz1:nz2),intent(inout) :: w
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
!-------------------------------------------------------------------------------
end subroutine extend_equal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function dist_point2plane(x0,y0,z0,A,B,C) result (L)
!-------------------------------------------------------------------------------
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
!-------------------------------------------------------------------------------
end function dist_point2plane
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine indx_locate_1d(x0,x,i1,i2,x1)
!-------------------------------------------------------------------------------
    real(SP),intent(in) :: x0
    real(SP),dimension(:),intent(in) :: x
    integer,intent(out) :: i1,i2
    real(SP),intent(out) :: x1
    integer :: i0
    !logical,intent(in) :: backward
      !integer :: indx(1)
      !indx=minloc(x,x>=x0)
      !i2=max(indx(1),2)
      !i1=i2-1
    
    i1=1; i2=size(x); x1=x0
    if (x0<=x(1)) then
       x1=x(1); i1=1; i2=i1+1
    elseif (x0>=x(i2)) then
       x1=x(i2); i1=i2-1
    else
      do while (i2-i1>1)
         i0=(i2-i1)/2+i1
         if (x(i0)==x0) then
            i1=i0;i2=i0+1; exit
         elseif (x(i0)>x0) then
            i2=i0;
         else
            i1=i0
         end if
      end do
    end if
!-------------------------------------------------------------------------------
end subroutine indx_locate_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine media_create(n_i,n_j,n_k)
!-------------------------------------------------------------------------------
    integer,intent(in) :: n_i,n_j,n_k
    character (len=SEIS_STRLEN) :: filenm
    integer,dimension(SEIS_GEO) ::subs,subc,subt
    filenm=media_fnm_get(n_i,n_j,n_k)
    call nfseis_grid3d_skel(filenm,nx,ny,nz,"media generated by seis3d_media" )
    subs=(/ nx1,ny1,nz1 /); subc=(/ nx,ny,nz /); subt=(/ 1,1,1 /)
    call nfseis_grid3d_attput(filenm,        &
         subs,subc,subt,                   &
         ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk, &
         nx1,nx2,ny1,ny2,nz1,nz2,nx,ny,nz, &
         ngi1,ngi2,ngj1,ngj2,ngk1,ngk2,    &
         ngx1,ngx2,ngy1,ngy2,ngz1,ngz2,    &
         (/ngi1,ngi2,ngj1,ngj2,ngk1,ngk2/) )
#ifdef WITHQS
    call nfseis_attput(filenm,'QsF0',QsF0)
    call nfseis_attput(filenm,'QsINF',QsINF)
    call nfseis_grid3d_addvar(filenm,'Qs')
#endif
    call nfseis_grid3d_addvar(filenm,'rho')
    call nfseis_grid3d_addvar(filenm,'C13')
    call nfseis_grid3d_addvar(filenm,'C66')
#if defined AnisVTI || defined AnisGene
    call nfseis_grid3d_addvar(filenm,'C11')
    call nfseis_grid3d_addvar(filenm,'C33')
    call nfseis_grid3d_addvar(filenm,'C44')
#endif
#ifdef AnisGene
    !call nfseis_grid3d_addvar(filenm,'C11')
     call nfseis_grid3d_addvar(filenm,'C12')
    !call nfseis_grid3d_addvar(filenm,'C13')
     call nfseis_grid3d_addvar(filenm,'C14')
     call nfseis_grid3d_addvar(filenm,'C15')
     call nfseis_grid3d_addvar(filenm,'C16')
     call nfseis_grid3d_addvar(filenm,'C22')
     call nfseis_grid3d_addvar(filenm,'C23')
     call nfseis_grid3d_addvar(filenm,'C24')
     call nfseis_grid3d_addvar(filenm,'C25')
     call nfseis_grid3d_addvar(filenm,'C26')
    !call nfseis_grid3d_addvar(filenm,'C33')
     call nfseis_grid3d_addvar(filenm,'C34')
     call nfseis_grid3d_addvar(filenm,'C35')
     call nfseis_grid3d_addvar(filenm,'C36')
    !call nfseis_grid3d_addvar(filenm,'C44')
     call nfseis_grid3d_addvar(filenm,'C45')
     call nfseis_grid3d_addvar(filenm,'C46')
     call nfseis_grid3d_addvar(filenm,'C55')
     call nfseis_grid3d_addvar(filenm,'C56')
#endif

    call nfseis_varput(filenm,'rho',rho,subs,subc,subt)
    call nfseis_varput(filenm,'C13',C13,subs,subc,subt)
    call nfseis_varput(filenm,'C66',C66,subs,subc,subt)
#if defined AnisVTI || defined AnisGene
    call nfseis_varput(filenm,'C11',C11,subs,subc,subt)
    call nfseis_varput(filenm,'C33',C33,subs,subc,subt)
    call nfseis_varput(filenm,'C44',C44,subs,subc,subt)
#endif
#ifdef AnisGene
    !call nfseis_varput(filenm,'C11',C11,subs,subc,subt)
     call nfseis_varput(filenm,'C12',C12,subs,subc,subt)
    !call nfseis_varput(filenm,'C13',C13,subs,subc,subt)
     call nfseis_varput(filenm,'C14',C14,subs,subc,subt)
     call nfseis_varput(filenm,'C15',C15,subs,subc,subt)
     call nfseis_varput(filenm,'C16',C16,subs,subc,subt)
     call nfseis_varput(filenm,'C22',C22,subs,subc,subt)
     call nfseis_varput(filenm,'C23',C23,subs,subc,subt)
     call nfseis_varput(filenm,'C24',C24,subs,subc,subt)
     call nfseis_varput(filenm,'C25',C25,subs,subc,subt)
     call nfseis_varput(filenm,'C26',C26,subs,subc,subt)
    !call nfseis_varput(filenm,'C33',C33,subs,subc,subt)
     call nfseis_varput(filenm,'C34',C34,subs,subc,subt)
     call nfseis_varput(filenm,'C35',C35,subs,subc,subt)
     call nfseis_varput(filenm,'C36',C36,subs,subc,subt)
    !call nfseis_varput(filenm,'C44',C44,subs,subc,subt)
     call nfseis_varput(filenm,'C45',C45,subs,subc,subt)
     call nfseis_varput(filenm,'C46',C46,subs,subc,subt)
     call nfseis_varput(filenm,'C55',C55,subs,subc,subt)
     call nfseis_varput(filenm,'C56',C56,subs,subc,subt)
#endif

#ifdef WITHQS
    call nfseis_varput(filenm,'Qs',Qs,subs,subc,subt)
#endif
!-------------------------------------------------------------------------------
end subroutine media_create
!-------------------------------------------------------------------------------

#ifdef MediaMPI
!-------------------------------------------------------------------------------
subroutine media_exchange
!-------------------------------------------------------------------------------
    integer :: ierr
    integer,dimension(MPI_STATUS_SIZE) :: istatus
!-- exchange
call MPI_SENDRECV(rho(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),111,        &
                  rho(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),111,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),112, &
                  rho(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),112, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),121,        &
                  rho(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),121,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),122, &
                  rho(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),122, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),131,        &
                  rho(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),131,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(rho(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),132, &
                  rho(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),132, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C66(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),211,        &
                  C66(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),211,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C66(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),212, &
                  C66(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),212, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C66(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),221,        &
                  C66(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),221,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C66(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),222, &
                  C66(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),222, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C66(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),231,        &
                  C66(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),231,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C66(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),232, &
                  C66(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),232, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C13(ni1  ,ny1,nz1),1,DTypeXL,neigid(1,1),311,        &
                  C13(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),311,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C13(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),312, &
                  C13(ni1-LenFDL  ,ny1,nz1),1,DTypeXL,neigid(1,1),312, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C13(nx1,nj1  ,nz1),1,DTypeYL,neigid(2,1),321,        &
                  C13(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),321,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C13(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),322, &
                  C13(nx1,nj1-LenFDL  ,nz1),1,DTypeYL,neigid(2,1),322, &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C13(nx1,ny1,nk1  ),1,DTypeZL,neigid(3,1),331,        &
                  C13(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),331,        &
                  SWMPI_COMM,istatus,ierr)
call MPI_SENDRECV(C13(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),332, &
                  C13(nx1,ny1,nk1-LenFDL  ),1,DTypeZL,neigid(3,1),332, &
                  SWMPI_COMM,istatus,ierr)
!-------------------------------------------------------------------------------
end subroutine media_exchange
!-------------------------------------------------------------------------------
#else
!-------------------------------------------------------------------------------
subroutine media_exchange
!-------------------------------------------------------------------------------
    character (len=SEIS_STRLEN) :: filenm
    integer,dimension(SEIS_GEO) :: subt
    integer,dimension(SEIS_GEO) ::         &
         subs_x1,subc_x1, subs_x2,subc_x2, &
         subs_y1,subc_y1, subs_y2,subc_y2, &
         subs_z1,subc_z1, subs_z2,subc_z2
    integer,dimension(LenFD) ::            &
         indx_x1,indx_x2,indx_y1,indx_y2,indx_z1,indx_z2
    integer :: n_i,n_j,n_k
    integer :: i,j,k
    
    subt=(/ 1,1,1 /)
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

do n_i=0,dims(1)-1
do n_j=0,dims(2)-1
do n_k=0,dims(3)-1
   write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
   call swmpi_change_fnm(n_i,n_j,n_k)
   call swmpi_set_gindx(n_i,n_j,n_k)
   call media_import(n_i,n_j,n_k)

   ! to x1
if (n_i>0) then
   call swmpi_change_fnm(n_i-1,n_j,n_k)
   filenm=media_fnm_get(n_i-1,n_j,n_k)
   call nfseis_varput(filenm,'rho',rho(indx_x1,:,:),subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'C66' ,C66(indx_x1,:,:), subs_x1,subc_x1,subt)
   call nfseis_varput(filenm,'C13',C13(indx_x1,:,:),subs_x1,subc_x1,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(indx_x1,:,:),subs_x1,subc_x1,subt)
#endif
end if
! to x2
if (n_i<dims(1)-1) then
   call swmpi_change_fnm(n_i+1,n_j,n_k)
   filenm=media_fnm_get(n_i+1,n_j,n_k)
   call nfseis_varput(filenm,'rho',rho(indx_x2,:,:),subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'C66', C66(indx_x2,:,:),subs_x2,subc_x2,subt)
   call nfseis_varput(filenm,'C13',C13(indx_x2,:,:),subs_x2,subc_x2,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(indx_x2,:,:),subs_x2,subc_x2,subt)
#endif
end if
! to y1
if (n_j>0) then
   call swmpi_change_fnm(n_i,n_j-1,n_k)
   filenm=media_fnm_get(n_i,n_j-1,n_k)
   call nfseis_varput(filenm,'rho',rho(:,indx_y1,:),subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'C66', C66(:,indx_y1,:), subs_y1,subc_y1,subt)
   call nfseis_varput(filenm,'C13',C13(:,indx_y1,:),subs_y1,subc_y1,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,indx_y1,:),subs_y1,subc_y1,subt)
#endif
end if
! to y2
if (n_j<dims(2)-1) then
   call swmpi_change_fnm(n_i,n_j+1,n_k)
   filenm=media_fnm_get(n_i,n_j+1,n_k)
   call nfseis_varput(filenm,'rho',rho(:,indx_y2,:),subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'C66', C66(:,indx_y2,:), subs_y2,subc_y2,subt)
   call nfseis_varput(filenm,'C13',C13(:,indx_y2,:),subs_y2,subc_y2,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,indx_y2,:),subs_y2,subc_y2,subt)
#endif
end if
! to k1
if (n_k>0) then
   call swmpi_change_fnm(n_i,n_j,n_k-1)
   filenm=media_fnm_get(n_i,n_j,n_k-1)
   call nfseis_varput(filenm,'rho', rho(:,:,indx_z1),subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'C66', C66(:,:,indx_z1), subs_z1,subc_z1,subt)
   call nfseis_varput(filenm,'C13',C13(:,:,indx_z1),subs_z1,subc_z1,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,:,indx_z1),subs_z1,subc_z1,subt)
#endif
end if
! to k2
if (n_k<dims(3)-1) then
   call swmpi_change_fnm(n_i,n_j,n_k+1)
   filenm=media_fnm_get(n_i,n_j,n_k+1)
   call nfseis_varput(filenm,'rho', rho(:,:,indx_z2),subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'C66', C66(:,:,indx_z2),  subs_z2,subc_z2,subt)
   call nfseis_varput(filenm,'C13',C13(:,:,indx_z2),subs_z2,subc_z2,subt)
#ifdef WITHQS
   call nfseis_varput(filenm,'Qs',Qs(:,:,indx_z2),subs_z2,subc_z2,subt)
#endif
end if
end do
end do
end do
!-------------------------------------------------------------------------------
end subroutine media_exchange
!-------------------------------------------------------------------------------
#endif

!-------------------------------------------------------------------------------
subroutine media_stept_check(n_i,n_j,n_k)
!-------------------------------------------------------------------------------
    integer,intent(in) :: n_i,n_j,n_k
    real(SP) :: Vp,dtLe,dtlocal
    integer :: i,j,k

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
#ifdef AnisVTI
       Vp=max( sqrt(C11(i,j,k)/rho(i,j,k) ), &
               sqrt(C33(i,j,k)/rho(i,j,k) ) )
#else
       Vp=sqrt( (C13(i,j,k)+2.0*C66(i,j,k))/rho(i,j,k) )
#endif
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
    !write(*,"(a,3i5,a,f,a,3i5,a,2f)") "   dtmax in thread", n_i,n_j,n_k, " is", dtmax,  &
    !    " located on", dtindx," with Vp and dL", dtmaxVp,dtmaxL
    write(*,"(a,3i5,a,f,a,3i5,a,2f)") "   dtmax in thread", dtnode, " is", dtmax,  &
        " located on", dtindx," with Vp and dL", dtmaxVp,dtmaxL
!-------------------------------------------------------------------------------
end subroutine media_stept_check
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine alloc_local
!-------------------------------------------------------------------------------
    allocate(gx(nx1:nx2,ny1:ny2,nz1:nz2)); gx=0.0_SP
    allocate(gy(nx1:nx2,ny1:ny2,nz1:nz2)); gy=0.0_SP
    allocate(gz(nx1:nx2,ny1:ny2,nz1:nz2)); gz=0.0_SP
!-------------------------------------------------------------------------------
end subroutine alloc_local
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine dealloc_media
!-------------------------------------------------------------------------------
    if (L3D%yes) then
       deallocate(L3D%x)
       deallocate(L3D%y)
       deallocate(L3D%h)
       deallocate(L3D%Vp)
       deallocate(L3D%Vs)
       deallocate(L3D%Dp)
    end if
!-------------------------------------------------------------------------------
end subroutine dealloc_media
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function anis_name2id(anisnm) result(id)
!-------------------------------------------------------------------------------
    character (len=*),intent(in) :: anisnm
    integer :: id
    select case (trim(anisnm))
    case ('ISO')
         id=SIG_ANIS_ISO
    case ('VTI')
         id=SIG_ANIS_VTI
    case ('TTI')
         id=SIG_ANIS_TTI
    case ('TRICLINIC')
         id=SIG_ANIS_TRICLINIC
    case default
         call error_except("out of anisotrpy type range "//trim(anisnm))
    end select
!-------------------------------------------------------------------------------
end function anis_name2id
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
function stif_name2id(stifnm) result(id)
!-------------------------------------------------------------------------------
    character (len=*),intent(in) :: stifnm
    integer :: id
    select case (trim(stifnm))
    case ('velocity')
         id=SIG_STIF_VEL
    case ('tensor')
         id=SIG_STIF_TEN
    case default
         call error_except("out of stiffness type range "//trim(stifnm))
    end select
!-------------------------------------------------------------------------------
end function stif_name2id
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine error_except(msg)
!-------------------------------------------------------------------------------
    character (len=*),intent(in) :: msg
#ifdef MediaMPI
    integer :: ierr
#endif
    print *, trim(msg)
#ifdef MediaMPI
    call MPI_ABORT(SWMPI_COMM,1,ierr)
#else
    stop 1
#endif
!-------------------------------------------------------------------------------
end subroutine error_except
!-------------------------------------------------------------------------------

end program seis3d_media

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:

