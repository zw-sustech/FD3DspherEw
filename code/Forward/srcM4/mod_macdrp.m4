divert(-1)dnl
dnl macro to generate mod_macdrp.F90
changequote([,])dnl
divert(0)dnl

module macdrp_mod

!-------------------------------------------------------------------------------
! Description:
!-------------------------------------------------------------------------------
!
! This module contains the variables and subroutines
! used in the DRP/opt MacCormack fd operator
!
!    Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
!    Copyright (C) 2006 Wei ZHANG
!
!-------------------------------------------------------------------------------
! Time stamp, log, version:
!-------------------------------------------------------------------------------
!
! $Date: 2010-03-02 19:54:21 -0500 (Tue, 02 Mar 2010) $
! $Revision: 94 $
! $LastChangedBy: zhangw $
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Include file:
!-------------------------------------------------------------------------------
#include "mod_macdrp.h"
!#define WATER

!-------------------------------------------------------------------------------
! Dependent modules:
!-------------------------------------------------------------------------------
    use constants_mod, only : SEIS_GEO
    use math_mod
    use para_mod
    use mpi
    use mpi_mod
    use media_mod
    use grid_mod
    use io_mod

!-------------------------------------------------------------------------------
! Public interface and variables:
!-------------------------------------------------------------------------------
    implicit none

    private

    !--- variables ---
    public ::                                      &
        Txx, Tyy, Txy, Vx, Vy, Tzz, Txz, Tyz, Vz,  &
        hTxx,hTyy,hTxy,hVx,hVy,hTzz,hTxz,hTyz,hVz

    !-- subroutine and function
    public ::                                      &
        macdrp_init,                               &
        macdrp_syn,                                &
        macdrp_mesg_init,                          &
        macdrp_destroy,                            &
        macdrp_LxF_LyF_LzF,                        &
        macdrp_LxB_LyB_LzB,                        &
        macdrp_LxB_LyB_LzF,                        &
        macdrp_LxF_LyF_LzB,                        &
        macdrp_LxB_LyF_LzF,                        &
        macdrp_LxF_LyB_LzB,                        &
        macdrp_LxF_LyB_LzF,                        &
        macdrp_LxB_LyF_LzB,                        &
        macdrp_RK_beg,                             &
        macdrp_RK_inn,                             &
        macdrp_RK_fin,                             &
        macdrp_check,                              &
        atten_graves,                              &
        macdrp_print_info,                         &
        macdrp_get_fnm_coupling

!-------------------------------------------------------------------------------
! Local interface
!-------------------------------------------------------------------------------
    interface macdrp_LxF_LyF_LzF
      module procedure in_LxF_LyF_LzF
    end interface
    interface macdrp_LxB_LyB_LzB
      module procedure in_LxB_LyB_LzB
    end interface
    interface macdrp_LxF_LyF_LzB
      module procedure in_LxF_LyF_LzB
    end interface
    interface macdrp_LxB_LyB_LzF
      module procedure in_LxB_LyB_LzF
    end interface
    interface macdrp_LxB_LyF_LzF
      module procedure in_LxB_LyF_LzF
    end interface
    interface macdrp_LxF_LyB_LzB
      module procedure in_LxF_LyB_LzB
    end interface
    interface macdrp_LxF_LyB_LzF
      module procedure in_LxF_LyB_LzF
    end interface
    interface macdrp_LxB_LyF_LzB
      module procedure in_LxB_LyF_LzB
    end interface

!-------------------------------------------------------------------------------
! Constant parameters in the include file:
!-------------------------------------------------------------------------------
DEFFDWET
DEFFDWET24
DEFFDWET22
DEFLDDRK2A
DEFLDDRK2B
DEFLDDRK4A
DEFLDDRK4B
HOCWETL
HOCWETR

!-------------------------------------------------------------------------------
! Module variables:
!-------------------------------------------------------------------------------
#ifdef MPIBuffered
    real(SP),dimension(:),allocatable ::        &
         BufX1,BufX2,BufY1,BufY2,BufZ1,BufZ2,   &
         RevX1,RevX2,RevY1,RevY2,RevZ1,RevZ2
    integer,parameter,private :: NREQ=4
    integer,private :: &
         NBufXL, NBufXS, &
         NBufYL, NBufYS, &
         NBufZL, NBufZS
#else
    !integer,parameter,private :: NREQ=36
    integer,parameter,private :: NREQ=24
#endif

    real(SP),dimension(:,:,:),allocatable ::        &
          Txx, Tyy, Txy, Vx, Vy, Tzz, Txz, Tyz, Vz, &
         hTxx,hTyy,hTxy,hVx,hVy,hTzz,hTxz,hTyz,hVz, &
         mTxx,mTyy,mTxy,mVx,mVy,mTzz,mTxz,mTyz,mVz, &
         tTxx,tTyy,tTxy,tVx,tVy,tTzz,tTxz,tTyz,tVz

#ifdef AnisGene
    real(SP),allocatable,public :: &
         matEh2Ev(:,:,:,:), matF2Ev(:,:,:,:)
#endif

    real(SP),dimension(:,:),allocatable,public :: &
         TxSrc,TySrc,TzSrc,                       &
         VxSrc,VySrc,VzSrc

!-- Far-Near coupling
    logical :: flag_couple_feature, flag_couple_node,flag_Far_read
    character (len=SEIS_STRLEN) :: pnm_coupling
    real(SP),dimension(:),allocatable ::    &
        Vx_Far_1D , Vy_Far_1D , Vz_Far_1D , &
        Txx_Far_1D, Tyy_Far_1D, Tzz_Far_1D, &
        Tyz_Far_1D, Txz_Far_1D, Txy_Far_1D
    integer,dimension(2*SEIS_GEO,2*SEIS_GEO) :: indx_Far_1D   !- index of first element in the 1D vector
    integer,dimension(2,SEIS_GEO) :: indx_couple_grid_center !- index of coupling grid without ghost index
    integer :: size_Far_1d     !- length of the 1D vector
    integer :: ntwin_Far_field !- time steps of the far field
    integer :: fid_Far         !- fid of the coupling input file
    integer,dimension(SEIS_GEO*2,SEIS_GEO*2+1),public :: scaindx !- index of scattering region
    logical,dimension(2,SEIS_GEO) :: flag_couple

!-- parameter of Runge-Kutta scheme
    real(SP),dimension(4),public :: firRKa,firRKb, secRKa,secRKb
    integer,dimension(SEIS_GEO*2,SEIS_GEO*2+1),public :: indx

!-- for MPI API
    integer :: ierr
    integer,dimension(MPI_STATUS_SIZE) :: istatus
    integer,dimension(NREQ) :: reqXB, reqXF, reqYB, reqYF, reqZB, reqZF
    integer,dimension(MPI_STATUS_SIZE,NREQ) :: reqstat

#ifdef VERBOSE
    integer fid_out
#endif

!===============================================================================
! subroutine and functions in this module
!===============================================================================

contains

!===============================================================================
subroutine macdrp_init(filenm)
!===============================================================================

!-------------------------------------------------------------------------------
!-- Description:
!-------------------------------------------------------------------------------
!--   Allocate module variables

!-------------------------------------------------------------------------------
!-- ToDo:
!-------------------------------------------------------------------------------
!--   Error handler

!-------------------------------------------------------------------------------
!-- Input arguments:
!-------------------------------------------------------------------------------
    character (len=*),intent(in) :: filenm

!-------------------------------------------------------------------------------
!-- Output arguments:
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!-- Local variables:
!-------------------------------------------------------------------------------
    integer :: ierr  !-- stat of allocate
    integer :: i,j,k,n

!-------------------------------------------------------------------------------
!-- Entry point: 
!-------------------------------------------------------------------------------
!-
    fid_Far = io_next_file_id()
    open(fid_Far,file=trim(filenm),status="old")
      call string_conf(fid_Far,1,'far_near_coupling',2,flag_couple_feature)
      if (flag_couple_feature) then
         call string_conf(fid_Far,1,'COUPLING_ROOT',2,pnm_coupling)
      end if
    close(fid_Far)
    if (flag_couple_feature) then
       open(fid_Far,file=macdrp_get_fnm_coupling(pnm_coupling,thisid(1),thisid(2),thisid(3)), &
                    form='unformatted', status='old')
       flag_Far_read=.false. !- the zero field assumed at zero time at local scale
    end if

    allocate( Txx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txx=0.0_SP
    allocate( Tyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tyy=0.0_SP
    allocate( Txy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txy=0.0_SP
    allocate( Vx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vx =0.0_SP
    allocate( Vy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vy =0.0_SP
    allocate( Tzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tzz=0.0_SP
    allocate( Txz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Txz=0.0_SP
    allocate( Tyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Tyz=0.0_SP
    allocate( Vz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr);  Vz =0.0_SP
    allocate(hTxx(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxx=0.0_SP
    allocate(hTyy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTyy=0.0_SP
    allocate(hTxy(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxy=0.0_SP
    allocate(hVx (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVx =0.0_SP
    allocate(hVy (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVy =0.0_SP
    allocate(hTzz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTzz=0.0_SP
    allocate(hTxz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTxz=0.0_SP
    allocate(hTyz(nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hTyz=0.0_SP
    allocate(hVz (nx1:nx2,ny1:ny2,nz1:nz2),stat=ierr); hVz =0.0_SP
    allocate(mTxx(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTxx=0.0_SP
    allocate(mTyy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTyy=0.0_SP
    allocate(mTxy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTxy=0.0_SP
    allocate(mVx (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mVx =0.0_SP
    allocate(mVy (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mVy =0.0_SP
    allocate(mTzz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTzz=0.0_SP
    allocate(mTxz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTxz=0.0_SP
    allocate(mTyz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mTyz=0.0_SP
    allocate(mVz (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); mVz =0.0_SP
    allocate(tTxx(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTxx=0.0_SP
    allocate(tTyy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTyy=0.0_SP
    allocate(tTxy(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTxy=0.0_SP
    allocate(tVx (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tVx =0.0_SP
    allocate(tVy (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tVy =0.0_SP
    allocate(tTzz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTzz=0.0_SP
    allocate(tTxz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTxz=0.0_SP
    allocate(tTyz(ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tTyz=0.0_SP
    allocate(tVz (ni1:ni2,nj1:nj2,nk1:nk2),stat=ierr); tVz =0.0_SP
    if (ierr>0) then
       print *, "can't allocate variable in macdrp_init"
       stop 1
    end if
    allocate(TxSrc(nx1:nx2,ny1:ny2),stat=ierr); TxSrc=0.0_SP
    allocate(TySrc(nx1:nx2,ny1:ny2),stat=ierr); TySrc=0.0_SP
    allocate(TzSrc(nx1:nx2,ny1:ny2),stat=ierr); TzSrc=0.0_SP
    allocate(VxSrc(nx1:nx2,ny1:ny2),stat=ierr); VxSrc=0.0_SP
    allocate(VySrc(nx1:nx2,ny1:ny2),stat=ierr); VySrc=0.0_SP
    allocate(VzSrc(nx1:nx2,ny1:ny2),stat=ierr); VzSrc=0.0_SP
#ifdef MPIBuffered
    allocate(BufX1(LenFD*nj*nk*6),stat=ierr); BufX1=0.0_SP
    allocate(BufX2(LenFD*nj*nk*6),stat=ierr); BufX2=0.0_SP
    allocate(BufY1(LenFD*ni*nk*6),stat=ierr); BufY1=0.0_SP
    allocate(BufY2(LenFD*ni*nk*6),stat=ierr); BufY2=0.0_SP
    allocate(BufZ1(LenFD*ni*nj*6),stat=ierr); BufZ1=0.0_SP
    allocate(BufZ2(LenFD*ni*nj*6),stat=ierr); BufZ2=0.0_SP
    allocate(RevX1(LenFD*nj*nk*6),stat=ierr); RevX1=0.0_SP
    allocate(RevX2(LenFD*nj*nk*6),stat=ierr); RevX2=0.0_SP
    allocate(RevY1(LenFD*ni*nk*6),stat=ierr); RevY1=0.0_SP
    allocate(RevY2(LenFD*ni*nk*6),stat=ierr); RevY2=0.0_SP
    allocate(RevZ1(LenFD*ni*nj*6),stat=ierr); RevZ1=0.0_SP
    allocate(RevZ2(LenFD*ni*nj*6),stat=ierr); RevZ2=0.0_SP
    NBufXL=nj*nk*LenFDL*6
    NBufXS=nj*nk*LenFDS*6
    NBufYL=ni*nk*LenFDL*6
    NBufYS=ni*nk*LenFDS*6
    NBufZL=ni*nj*LenFDL*6
    NBufZS=ni*nj*LenFDS*6
#endif

!-- inner region without points affected by MPI layers
    indx(:,SEIS_GEO*2+1)=(/ ni1+LenFD,ni2-LenFD, &
                            nj1+LenFD,nj2-LenFD, &
                            nk1+LenFD,nk2-LenFD /)
!-- boundary layer calculated after message exchange
    indx(:,SEIS_GEO*2  )=(/ ni1,ni2,nj1,nj2,nk2-LenFD+1,nk2 /) ! z2
    indx(:,SEIS_GEO*2-1)=(/ ni1,ni2,nj1,nj2,nk1,nk1+LenFD-1 /) ! z1
    indx(:,3)=(/ ni1,ni2,nj1,nj1+LenFD-1,nk1+LenFD,nk2-LenFD /) ! y1
    indx(:,4)=(/ ni1,ni2,nj2-LenFD+1,nj2,nk1+LenFD,nk2-LenFD /) ! y2
    indx(:,1)=(/ ni1,ni1+LenFD-1,nj1+LenFD,nj2-LenFD,nk1+LenFD,nk2-LenFD /) ! x1
    indx(:,2)=(/ ni2-LenFD+1,ni2,nj1+LenFD,nj2-LenFD,nk1+LenFD,nk2-LenFD /) ! x2
!-- initialize scattering region to non exist
    scaindx(1:5:2,:)=0; scaindx(2:6:2,:)=-1

#ifdef VERBOSE
!-- print index before coupling
    write(myid+10000,*) "indx before coupling"
    write(myid+10000,*) indx(:,1)
    write(myid+10000,*) indx(:,2)
    write(myid+10000,*) indx(:,3)
    write(myid+10000,*) indx(:,4)
    write(myid+10000,*) indx(:,5)
    write(myid+10000,*) indx(:,6)
    write(myid+10000,*) indx(:,7)
#endif

!-- rk coefficient
    firRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
    firRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
    !secRKa=(/ RK4a2, RK4a3, RK4a4, 0.0 /)
    !secRKb=(/ RK4b1, RK4b2, RK4b3, RK4b4 /)
    secRKa=(/ RK2a2, 0.0, 0.0, 0.0 /)
    secRKb=(/ RK2b1, RK2b2, 0.0, 0.0 /)

#ifdef AnisGene
!-- mat to convert V,z
    call coef_Eh2Ev
#endif

!-- Far-Near coupling
    flag_couple=.false.

    if (flag_couple_feature) then
    read(fid_Far) flag_couple_node,size_Far_1D,indx_Far_1D,ntwin_Far_field
    if ( flag_couple_node ) then
    !-- check grid size
        read(fid_Far) i,j,k
        if ( i /= ni .or. j /= nj .or. k /= nk ) then
           write(*,*) "length in coupling file:",i,j,k
           write(*,*) "grid size of modeling:",ni,nj,nk
           call swmpi_except("Error: Inconsistent of length in coupling file and model")
        end if
    !-- initialize coupling center grid to -1 meaning no coupling
        indx_couple_grid_center=-1
    !-- shift indx_Far_1D due to ghost points
        indx_Far_1D=indx_Far_1D+LenFD
    !-- calculate coupling center grid
        if (indx_Far_1D(2,1)>indx_Far_1D(1,1)) then
              indx_couple_grid_center(1,1)=(indx_Far_1D(1,1)+indx_Far_1D(2,1)+1)/2
              flag_couple(1,1)=.true.
        else
              indx_couple_grid_center(1,1)=ni1
        end if
        if (indx_Far_1D(2,2)>indx_Far_1D(1,2)) then
              indx_couple_grid_center(2,1)=(indx_Far_1D(1,2)+indx_Far_1D(2,2)  )/2
              flag_couple(2,1)=.true.
        else
              indx_couple_grid_center(2,1)=nj1
        end if
        if (indx_Far_1D(4,3)>indx_Far_1D(3,3)) then
              indx_couple_grid_center(1,2)=(indx_Far_1D(3,3)+indx_Far_1D(4,3)+1)/2
              flag_couple(1,2)=.true.
        else
              indx_couple_grid_center(1,2)=nj1
        end if
        if (indx_Far_1D(4,4)>indx_Far_1D(3,4)) then
              indx_couple_grid_center(2,2)=(indx_Far_1D(3,4)+indx_Far_1D(4,4)  )/2
              flag_couple(2,2)=.true.
        else
              indx_couple_grid_center(2,2)=nj2
        end if
        if (indx_Far_1D(6,5)>indx_Far_1D(5,5)) then
              indx_couple_grid_center(1,3)=(indx_Far_1D(5,5)+indx_Far_1D(6,5)+1)/2
              flag_couple(1,3)=.true.
        else
              indx_couple_grid_center(1,3)=nk1
        end if
        if (indx_Far_1D(6,6)>indx_Far_1D(5,6)) then
              indx_couple_grid_center(2,3)=(indx_Far_1D(5,6)+indx_Far_1D(6,6)  )/2
              flag_couple(2,3)=.true.
        else
              indx_couple_grid_center(2,3)=nk2
        end if

    !-- allocate
        allocate( Vx_Far_1D(size_Far_1D),stat=ierr);  Vx_Far_1D = 0.0_SP
        allocate( Vy_Far_1D(size_Far_1D),stat=ierr);  Vy_Far_1D = 0.0_SP
        allocate( Vz_Far_1D(size_Far_1D),stat=ierr);  Vz_Far_1D = 0.0_SP
        allocate(Txx_Far_1D(size_Far_1D),stat=ierr); Txx_Far_1D = 0.0_SP
        allocate(Tyy_Far_1D(size_Far_1D),stat=ierr); Tyy_Far_1D = 0.0_SP
        allocate(Tzz_Far_1D(size_Far_1D),stat=ierr); Tzz_Far_1D = 0.0_SP
        allocate(Tyz_Far_1D(size_Far_1D),stat=ierr); Tyz_Far_1D = 0.0_SP
        allocate(Txz_Far_1D(size_Far_1D),stat=ierr); Txz_Far_1D = 0.0_SP
        allocate(Txy_Far_1D(size_Far_1D),stat=ierr); Txy_Far_1D = 0.0_SP
        !if ( any(indx_couple_grid_center-LenFD/=0 .and. indx_couple_grid_center-LenFD<=0)) then
        !   call swmpi_except("coupling index should > LenFD")
        !end if

    !-- z1 is scattering margin
        if (flag_couple(1,3)) then
            scaindx(5,5)=nk1;
            scaindx(6,5)=indx_couple_grid_center(1,3)-1
            indx   (5,5)=scaindx(6,5)+1;
            indx   (6,5)=indx(5,5)+LenFD-1
            scaindx(1:4,5)=(/ ni1,ni2,nj1,nj2 /)  !- cover entire
            scaindx(5,1:4)=scaindx(6,5)+1 !- hori scat above lower scat
            indx   (5,1:4)=indx(6,5)+1    !- horizonal MPI layer above lower-MPI
        else
            scaindx(5,5)= 0;
            scaindx(6,5)=-1
            indx   (5,5)=nk1;
            indx   (6,5)=indx(5,5)+LenFD-1
            scaindx(5,1:4)=nk1       !- hori scat starts from bottom
            indx (5,1:4)=indx(6,5)+1 !- hori MPI layer above lower-MPI
        end if
        indx(5,7)=indx(6,5)+1 !- main always above lower-MPI layers
    !-- z2 is scattering margin
        if (flag_couple(2,3)) then
            scaindx(5,6)=indx_couple_grid_center(2,3)+1
            scaindx(6,6)=nk2;
            indx   (6,6)=scaindx(5,6)-1;
            indx   (5,6)=indx(6,6)-LenFD+1
            scaindx(1:4,6)=(/ ni1,ni2,nj1,nj2 /) !- cover entire
            scaindx(6,1:4)=scaindx(5,6)-1        !- hori scat below upper scat
            indx   (6,1:4)=indx(5,6)-1           !- hori MPI layer below upper-MPI
        else
            scaindx(5,6)= 0;
            scaindx(6,6)=-1
            indx   (6,6)=nk2;
            indx   (5,6)=indx(6,6)-LenFD+1
            scaindx(6,1:4)=nk2                   !- hori scat reaches top
            indx (6,1:4)=indx(5,6)-1             !- hori MPI layer below upper-MPI
        end if
        indx(6,7)=indx(5,6)-1 !- main always below upper-MPI layers

    !-- in y1 y2, indx(3:4,5:7) should be set
    !-- y1 is scattering
        if (flag_couple(1,2)) then
            scaindx(3,3)=nj1
            scaindx(4,3)=indx_couple_grid_center(1,2)-1
            indx(3,3)=scaindx(4,3)+1
            indx(4,3)=indx(3,3)+LenFD-1
            scaindx(1:2,3)=(/ ni1,ni2 /)  !- y scat cover overlap with x scat
            scaindx(3,1:2)=scaindx(4,3)+1 !- x1,x2 scatter region exclude y1
            indx(3,1:2)=indx(4,3)+1       !- x MPI layer exclude y
            indx(3,5)=scaindx(4,3)+1      !- vert MPI layer right to y1 scat
            indx(3,6)=scaindx(4,3)+1      !- vert MPI layer right to y1 scat
        else
            scaindx(3,3)=0
            scaindx(4,3)=-1
            indx(3,3)=nj1
            indx(4,3)=nj1+LenFD
            scaindx(3,1:2)=nj1        !- x1,x2 scatter region starts from y1
            indx(3,1:2)=indx(4,3)+1   !- x MPI layer exclude y
            indx(3,5)=nj1             !- vert MPI layer starts y1
            indx(3,6)=nj1             !- vert MPI layer starts y1
        end if
        indx(3,7)=indx(4,3)+1 !- main always right to left-MPI layers
    !-- y2 is scattering
        if (flag_couple(2,2)) then
            scaindx(3,4)=indx_couple_grid_center(2,2)+1
            scaindx(4,4)=nj2
            indx(4,4)=scaindx(3,4)-1
            indx(3,4)=indx(4,4)-LenFD+1
            scaindx(1:2,4)=(/ ni1,ni2 /)  !- y scat cover overlap with x scat
            scaindx(4,1:2)=scaindx(3,4)-1 !- x1,x2 scatter region exclude y2
            indx(4,1:2)=indx(3,4)-1       ! x MPI layer exclude y
            indx(4,5)=scaindx(3,4)-1      !- vert MPI layer left to y2 scat
            indx(4,6)=scaindx(3,4)-1      !- vert MPI layer left to y2 scat
        else
            scaindx(3,4)=0
            scaindx(4,4)=-1
            indx(3,4)=nj2-LenFD+1
            indx(4,4)=nj2
            scaindx(4,1:2)=nj2  !- x1,x2 scatter region reaches y2
            indx(4,1:2)=indx(3,4)-1  ! x MPI layer exclude y
            indx(4,5)=nj2      !- vert MPI layer reaches y2
            indx(4,6)=nj2      !- vert MPI layer reaches y2
        end if
        indx(4,7)=indx(3,4)-1 !- main always left to right-MPI layers

    !-- in x1 x2, indx(1:2,3:7) should be set
    !-- x1 is scattering
        if (flag_couple(1,1)) then
            scaindx(1,1)=ni1
            scaindx(2,1)=indx_couple_grid_center(1,1)-1
            indx(1,1)=scaindx(2,1)+1
            indx(2,1)=scaindx(2,1)+LenFD
            indx(1,3)=scaindx(2,1)+1
            indx(1,4)=scaindx(2,1)+1
            indx(1,5)=scaindx(2,1)+1
            indx(1,6)=scaindx(2,1)+1
        else
            scaindx(1,1)=0
            scaindx(2,1)=-1
            indx(1,1)=ni1
            indx(2,1)=ni1+LenFD
            indx(1,3)=ni1
            indx(1,4)=ni1
            indx(1,5)=ni1
            indx(1,6)=ni1
        end if
        indx(1,7)=indx(2,1)+1 !- main always right to left-MPI layers
    !-- x2 is scattering
        if (flag_couple(2,1)) then
            scaindx(1,2)=indx_couple_grid_center(2,1)+1
            scaindx(2,2)=ni2
            indx(2,2)=scaindx(1,2)-1
            indx(1,2)=indx(2,2)-LenFD+1
            indx(2,3)=scaindx(1,2)-1
            indx(2,4)=scaindx(1,2)-1
            indx(2,5)=scaindx(1,2)-1
            indx(2,6)=scaindx(1,2)-1
        else
            scaindx(1,2)=0
            scaindx(2,2)=-1
            indx(1,2)=ni2-LenFD+1
            indx(2,2)=ni2
            indx(2,3)=ni2
            indx(2,4)=ni2
            indx(2,5)=ni2
            indx(2,6)=ni2
        end if
        indx(2,7)=indx(1,2)-1 !- main always left to right-MPI layers

    !!-- main inner domain excludes MPI buffer and scattering margin
    !    do n = 1 , SEIS_GEO
    !    !-- negative face
    !        if (indx_couple_grid_center(1,n)>0) then
    !           indx(2*(n-1)+1,SEIS_GEO*2+1) = indx_couple_grid_center(1,n)+LenFD+LenFD !- ghost+MPI
    !        end if
    !    !-- positive face
    !        if (indx_couple_grid_center(2,n)>0) then
    !           indx(2*(n-1)+2,SEIS_GEO*2+1) = indx_couple_grid_center(2,n)+LenFD-LenFD !- ghost-MPI
    !        end if
    !    end do
    !!-- main MPI buffer surrounds main inner region and connects with scattering region
    !    indx(:,6)=(/ indx(1,7)-LenFD,indx(2,7)+LenFD,indx(3,7)-LenFD,indx(4,7)+LenFD,indx(6,7)+1,indx(6,7)+LenFD /) ! z2
    !    indx(:,5)=(/ indx(1,7)-LenFD,indx(2,7)+LenFD,indx(3,7)-LenFD,indx(4,7)+LenFD,indx(5,7)-LenFD,indx(5,7)-1 /) ! z1
    !    indx(:,4)=(/ indx(1,7)-LenFD,indx(2,7)+LenFD,indx(4,7)+1,indx(4,7)+LenFD,indx(5,7),indx(6,7) /) ! y2
    !    indx(:,3)=(/ indx(1,7)-LenFD,indx(2,7)+LenFD,indx(3,7)-LenFD,indx(3,7)-1,indx(5,7),indx(6,7) /) ! y1
    !    indx(:,2)=(/ indx(2,7)+1,indx(2,7)+LenFD,indx(3,7),indx(4,7),indx(5,7),indx(6,7) /) ! x2
    !    indx(:,1)=(/ indx(1,7)-LenFD,indx(1,7)-1,indx(3,7),indx(4,7),indx(5,7),indx(6,7) /) ! x1
    !!-- scattering region index
    !    scaindx(:,6)=(/ ni1,ni2,nj1,nj2,indx(6,6)+1,nk2 /) ! z2
    !    scaindx(:,5)=(/ ni1,ni2,nj1,nj2,nk1,indx(5,5)-1 /) ! z1
    !    scaindx(:,4)=(/ ni1,ni2,indx(4,4)+1,nj2,indx(5,5),indx(6,6) /) ! y2
    !    scaindx(:,3)=(/ ni1,ni2,nj1,indx(3,3)-1,indx(5,5),indx(6,6) /) ! y1
    !    scaindx(:,2)=(/ indx(2,2)+1,ni2,indx(3,3),indx(4,4),indx(5,5),indx(6,6) /) ! x2
    !    scaindx(:,1)=(/ ni1,indx(1,1)-1,indx(3,3),indx(4,4),indx(5,5),indx(6,6) /) ! x1
    end if
    end if

#ifdef VERBOSE
!-- print index after coupling
    write(myid+10000,*) "flag_couple_feature=",flag_couple_feature
    write(myid+10000,*) "flag_couple_node=",flag_couple_node
    write(myid+10000,*) "size_Far_1D=",size_Far_1D
    write(myid+10000,*) "indx_Far_1D=",indx_Far_1D
    write(myid+10000,*) "indx_couple_grid_center=",indx_couple_grid_center
    write(myid+10000,*) "indx after coupling"
    write(myid+10000,*) indx(:,1)
    write(myid+10000,*) indx(:,2)
    write(myid+10000,*) indx(:,3)
    write(myid+10000,*) indx(:,4)
    write(myid+10000,*) indx(:,5)
    write(myid+10000,*) indx(:,6)
    write(myid+10000,*) indx(:,7)
    write(myid+10000,*) "scattering region index"
    write(myid+10000,*) scaindx(:,1)
    write(myid+10000,*) scaindx(:,2)
    write(myid+10000,*) scaindx(:,3)
    write(myid+10000,*) scaindx(:,4)
    write(myid+10000,*) scaindx(:,5)
    write(myid+10000,*) scaindx(:,6)
#endif

#ifdef VERBOSE
    fid_out=9050
    open(fid_out,                                                                      &
         file='log_maxval_'//trim(set_mpi_subfix(thisid(1),thisid(2),thisid(3)))//'.dat', &
         status='unknown')
#endif
!===============================================================================
end subroutine macdrp_init
!===============================================================================

#ifdef AnisGene
!===============================================================================
subroutine coef_Eh2Ev
!===============================================================================

!-------------------------------------------------------------------------------
!-- Description:
!-------------------------------------------------------------------------------
!--   Coefficent to convert strain at the surface

!-------------------------------------------------------------------------------
!-- ToDo:
!-------------------------------------------------------------------------------
!--

!-------------------------------------------------------------------------------
!-- Input arguments:
!-------------------------------------------------------------------------------
!--

!-------------------------------------------------------------------------------
!-- Output arguments:
!-------------------------------------------------------------------------------
!--

!-------------------------------------------------------------------------------
!-- Global variables or function need:
!-------------------------------------------------------------------------------
    use math_mod, only : invert

!-------------------------------------------------------------------------------
!-- Local variables:
!-------------------------------------------------------------------------------
    real(SP),dimension(SEIS_GEO,SEIS_GEO) :: A,B     !-- for invert
    integer :: i,j,k                                 !-- loop variable
    integer :: ierr                                  !-- stat of allocate

!-------------------------------------------------------------------------------
!-- Entry point: 
!-------------------------------------------------------------------------------
!-

  !---- only needed for free surface ----
  if ( .not. freenode ) return

  allocate(matEh2Ev(SEIS_GEO,SEIS_GEO,ni1:ni2,nj1:nj2),stat=ierr); matEh2Ev=0.0_SP
  allocate(matF2Ev(SEIS_GEO,SEIS_GEO,ni1:ni2,nj1:nj2),stat=ierr); matF2Ev=0.0_SP

  k=nk2
  do j=nj1,nj2
  do i=ni1,ni2
     !----
     A(1,:)=(/ C33(i,j,k),C34(i,j,k),C35(i,j,k) /)
     A(2,:)=(/ C34(i,j,k),C44(i,j,k),C45(i,j,k) /)
     A(3,:)=(/ C35(i,j,k),C45(i,j,k),C55(i,j,k) /)
     call invert(A)
     !----
     B(1,:)=(/ C13(i,j,k),C23(i,j,k),C36(i,j,k) /)
     B(2,:)=(/ C14(i,j,k),C24(i,j,k),C46(i,j,k) /)
     B(3,:)=(/ C15(i,j,k),C25(i,j,k),C56(i,j,k) /)
     B=-B
     !----
     matEh2Ev(1:SEIS_GEO,1:SEIS_GEO,i,j)=matmul(A,B)
     matF2Ev (1,1:SEIS_GEO,i,j) = (/ A(1,1),A(1,2),A(1,3) /)
     matF2Ev (2,1:SEIS_GEO,i,j) = (/ A(2,1),A(2,2),A(2,3) /)
     matF2Ev (3,1:SEIS_GEO,i,j) = (/ A(3,1),A(3,2),A(3,3) /)
  end do
  end do

!===============================================================================
end subroutine coef_Eh2Ev
!===============================================================================
#endif

!===============================================================================
subroutine macdrp_destroy  !-- deallocate module variable, not complete
!===============================================================================
    deallocate( Txx, Tyy, Txy, Vx, Vy)
    deallocate(hTxx,hTyy,hTxy,hVx,hVy)
    deallocate(mTxx,mTyy,mTxy,mVx,mVy)
    deallocate(tTxx,tTyy,tTxy,tVx,tVy)
!-- Far-Near coupling
    if ( allocated( Vx_Far_1D) ) deallocate( Vx_Far_1D)
    if ( allocated( Vy_Far_1D) ) deallocate( Vy_Far_1D)
    if ( allocated( Vz_Far_1D) ) deallocate( Vz_Far_1D)
    if ( allocated(Txx_Far_1D) ) deallocate(Txx_Far_1D)
    if ( allocated(Tyy_Far_1D) ) deallocate(Tyy_Far_1D)
    if ( allocated(Tzz_Far_1D) ) deallocate(Tzz_Far_1D)
    if ( allocated(Tyz_Far_1D) ) deallocate(Tyz_Far_1D)
    if ( allocated(Txz_Far_1D) ) deallocate(Txz_Far_1D)
    if ( allocated(Txy_Far_1D) ) deallocate(Txy_Far_1D)
    close(fid_Far)
#ifdef VERBOSE
    close(fid_out)
#endif
!===============================================================================
end subroutine macdrp_destroy
!===============================================================================

!===============================================================================
subroutine macdrp_check(ntime)  !-- check overflow at ntime step
!===============================================================================
    integer,intent(in) :: ntime
    real(SP) :: V1,V2,V3,T11,T22,T33,T12,T13,T23,W
    integer ierr
#ifndef CheckOverFlow
    return
#endif
    if (mod(ntime,1)==0) then
        V1=maxval(abs(Vx))
        V2=maxval(abs(Vy))
        V3=maxval(abs(Vz))
       T11=maxval(abs(Txx))
       T22=maxval(abs(Tyy))
       T33=maxval(abs(Tzz))
       T12=maxval(abs(Txy))
       T13=maxval(abs(Txz))
       T23=maxval(abs(Tyz))
#ifdef VERBOSE
       write(fid_out,'(i5,9es12.5)') ntime, V1,V2,V3,T11,T22,T33,T12,T13,T23
#endif
       W=max(V1,V2,V3,T11,T22,T33,T12,T13,T23)
       if (W>=huge(1.0)) then
          print *, "Overflow error: "
          write(*,"(i5,i3.2,2(i2.2),9(es12.5,3i5))") ntime,thisid(1),thisid(2),thisid(3), &
             V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
            T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
            T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))
#ifdef VERBOSE
          write(fid_out,"(i5,9(es12.5,3i5))") ntime,                                      &
             V1, maxloc(abs(Vx)), V2, maxloc(abs(Vy)), V3, maxloc(abs(Vz)),               &
            T11, maxloc(abs(Txx)), T22, maxloc(abs(Tyy)), T33, maxloc(abs(Tzz)),          &
            T12, maxloc(abs(Txy)), T13, maxloc(abs(Txz)), T23, maxloc(abs(Tyz))
#endif
          call MPI_ABORT(SWMPI_COMM,1,ierr)
       end if
    end if
!===============================================================================
end subroutine macdrp_check
!===============================================================================

!===============================================================================
subroutine macdrp_syn
!===============================================================================
    integer i,j,k
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
     mTxx(i,j,k)=Txx(i,j,k)
     mTyy(i,j,k)=Tyy(i,j,k)
     mTxy(i,j,k)=Txy(i,j,k)
     mVx (i,j,k)=Vx (i,j,k)
     mVy (i,j,k)=Vy (i,j,k)
     mTzz(i,j,k)=Tzz (i,j,k)
     mTxz(i,j,k)=Txz (i,j,k)
     mTyz(i,j,k)=Tyz (i,j,k)
     mVz (i,j,k)=Vz  (i,j,k)
    end do
    end do
    end do
!$OMP END PARALLEL DO
!===============================================================================
end subroutine macdrp_syn
!===============================================================================

!===============================================================================
subroutine macdrp_RK_beg(rka,rkb)  !--- first step of RK scheme
!===============================================================================
    real(SP),intent(in) :: rka,rkb
    real(SP) :: a,b
    integer i,j,k
    a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
       Txx(i,j,k)=mTxx(i,j,k)+a*hTxx(i,j,k)
       Tyy(i,j,k)=mTyy(i,j,k)+a*hTyy(i,j,k)
       Tzz(i,j,k)=mTzz(i,j,k)+a*hTzz(i,j,k)
       Txy(i,j,k)=mTxy(i,j,k)+a*hTxy(i,j,k)
       Txz(i,j,k)=mTxz(i,j,k)+a*hTxz(i,j,k)
       Tyz(i,j,k)=mTyz(i,j,k)+a*hTyz(i,j,k)
       Vx (i,j,k)=mVx (i,j,k)+a*hVx (i,j,k)
       Vy (i,j,k)=mVy (i,j,k)+a*hVy (i,j,k)
       Vz (i,j,k)=mVz (i,j,k)+a*hVz (i,j,k)
    
       tTxx(i,j,k)=mTxx(i,j,k)+b*hTxx(i,j,k)
       tTyy(i,j,k)=mTyy(i,j,k)+b*hTyy(i,j,k)
       tTzz(i,j,k)=mTzz(i,j,k)+b*hTzz(i,j,k)
       tTxy(i,j,k)=mTxy(i,j,k)+b*hTxy(i,j,k)
       tTxz(i,j,k)=mTxz(i,j,k)+b*hTxz(i,j,k)
       tTyz(i,j,k)=mTyz(i,j,k)+b*hTyz(i,j,k)
       tVx (i,j,k)=mVx (i,j,k)+b*hVx (i,j,k)
       tVy (i,j,k)=mVy (i,j,k)+b*hVy (i,j,k)
       tVz (i,j,k)=mVz (i,j,k)+b*hVz (i,j,k)
    end do
    end do
    end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
    call free_charac
    call free_extrap
#endif
!===============================================================================
end subroutine macdrp_RK_beg
!===============================================================================

!===============================================================================
subroutine macdrp_RK_inn(rka,rkb)  !-- inner step of RK scheme
!===============================================================================
    real(SP),intent(in) :: rka,rkb
    real(SP) :: a,b
    integer i,j,k
    a=rka*stept; b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
       Txx(i,j,k)=mTxx(i,j,k)+a*hTxx(i,j,k)
       Tyy(i,j,k)=mTyy(i,j,k)+a*hTyy(i,j,k)
       Tzz(i,j,k)=mTzz(i,j,k)+a*hTzz(i,j,k)
       Txy(i,j,k)=mTxy(i,j,k)+a*hTxy(i,j,k)
       Txz(i,j,k)=mTxz(i,j,k)+a*hTxz(i,j,k)
       Tyz(i,j,k)=mTyz(i,j,k)+a*hTyz(i,j,k)
       Vx (i,j,k)=mVx (i,j,k)+a*hVx (i,j,k)
       Vy (i,j,k)=mVy (i,j,k)+a*hVy (i,j,k)
       Vz (i,j,k)=mVz (i,j,k)+a*hVz (i,j,k)
    
       tTxx(i,j,k)=tTxx(i,j,k)+b*hTxx(i,j,k)
       tTyy(i,j,k)=tTyy(i,j,k)+b*hTyy(i,j,k)
       tTzz(i,j,k)=tTzz(i,j,k)+b*hTzz(i,j,k)
       tTxy(i,j,k)=tTxy(i,j,k)+b*hTxy(i,j,k)
       tTxz(i,j,k)=tTxz(i,j,k)+b*hTxz(i,j,k)
       tTyz(i,j,k)=tTyz(i,j,k)+b*hTyz(i,j,k)
       tVx (i,j,k)=tVx (i,j,k)+b*hVx (i,j,k)
       tVy (i,j,k)=tVy (i,j,k)+b*hVy (i,j,k)
       tVz (i,j,k)=tVz (i,j,k)+b*hVz (i,j,k)
    end do
    end do
    end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
    call free_charac
    call free_extrap
#endif
!===============================================================================
end subroutine macdrp_RK_inn
!===============================================================================

!===============================================================================
subroutine macdrp_RK_fin(rkb)  !-- final step of RK scheme
!===============================================================================
    real(SP),intent(in) :: rkb
    real(SP) :: b
    integer i,j,k
    b=rkb*stept
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
       Txx(i,j,k)=tTxx(i,j,k)+b*hTxx(i,j,k)
       Tyy(i,j,k)=tTyy(i,j,k)+b*hTyy(i,j,k)
       Tzz(i,j,k)=tTzz(i,j,k)+b*hTzz(i,j,k)
       Txy(i,j,k)=tTxy(i,j,k)+b*hTxy(i,j,k)
       Txz(i,j,k)=tTxz(i,j,k)+b*hTxz(i,j,k)
       Tyz(i,j,k)=tTyz(i,j,k)+b*hTyz(i,j,k)
       Vx (i,j,k)=tVx (i,j,k)+b*hVx (i,j,k)
       Vy (i,j,k)=tVy (i,j,k)+b*hVy (i,j,k)
       Vz (i,j,k)=tVz (i,j,k)+b*hVz (i,j,k)
    end do
    end do
    end do
!$OMP END PARALLEL DO
#ifdef CondFreeCharac
    call free_charac
    call free_extrap
#endif
!===============================================================================
end subroutine macdrp_RK_fin
!===============================================================================

!===============================================================================
subroutine atten_graves  !-- attenuation effect
!===============================================================================
    integer :: i,j,k
    real(SP) :: Qatt
!$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k)
    do k=nk1,nk2
    do j=nj1,nj2
    do i=ni1,ni2
    if (Qs(i,j,k)<QsINF) then
       ! Qs is Q value not attenuation coefficient now
       Qatt=exp((-PI*QsF0*stept)/Qs(i,j,k))
       Txx(i,j,k)=Txx(i,j,k)*Qatt
       Tyy(i,j,k)=Tyy(i,j,k)*Qatt
       Tzz(i,j,k)=Tzz(i,j,k)*Qatt
       Txy(i,j,k)=Txy(i,j,k)*Qatt
       Txz(i,j,k)=Txz(i,j,k)*Qatt
       Tyz(i,j,k)=Tyz(i,j,k)*Qatt
       Vx (i,j,k)=Vx (i,j,k)*Qatt
       Vy (i,j,k)=Vy (i,j,k)*Qatt
       Vz (i,j,k)=Vz (i,j,k)*Qatt
    end if
    end do
    end do
    end do
!$OMP END PARALLEL DO
!===============================================================================
end subroutine atten_graves
!===============================================================================

!===============================================================================
!--- 3D bias FD of the MacCormack scheme
!===============================================================================

divert([-1])
define([_WRAPLXYZ],
[dnl
!-------------------------------------------------------------------------------
subroutine in_Lx[$1]_Ly[$2]_Lz[$3]
!-------------------------------------------------------------------------------
    integer n
#ifdef MPIBARRIER
    integer ierr
    call MPI_BARRIER(SWMPI_COMM,ierr)
#endif
    !Txx=real(myid)
#ifdef CondFreeTIMG
    if (freenode) call free_timg
#endif
#ifdef CondFreeVEXT
    if (freenode) call free_vext
#endif
    n=SEIS_GEO*2+1
#ifdef MPIBuffered
!-- pack message
    call fill_buff_Lx[$1]
    call fill_buff_Ly[$2]
    call fill_buff_Lz[$3]
#endif
!-- starting message exchanging
    call MPI_STARTALL(NREQ,reqX[$1],ierr)
    call MPI_STARTALL(NREQ,reqY[$2],ierr)
    call MPI_STARTALL(NREQ,reqZ[$3],ierr)
!-- read far field every two RK stages
    if (flag_Far_read) then
       call read_farfield
       flag_Far_read=.false.
    else
       flag_Far_read=.true.
    end if
!-- add far field in scatter region except MPI layers
    call add_far2scat_main
!-- calculation on main region
    call Lx[$1]_Ly[$2]_Lz[$3]( indx(1,n),indx(2,n), &
                      indx(3,n),indx(4,n), &
                      indx(5,n),indx(6,n) )
!-- wait the finish of message
    call MPI_WAITALL(NREQ,reqX[$1],reqstat,ierr)
    call MPI_WAITALL(NREQ,reqY[$2],reqstat,ierr)
    call MPI_WAITALL(NREQ,reqZ[$3],reqstat,ierr)
#ifdef MPIBuffered
!-- unpack message
    call recv_buff_Lx[$1]
    call recv_buff_Ly[$2]
    call recv_buff_Lz[$3]
#endif
!-- add far field in scatter region in MPI layers
    call add_far2scat_boundary
!-- calculation on boundary region inside coupling
    do n=1,SEIS_GEO*2
    call Lx[$1]_Ly[$2]_Lz[$3]( indx(1,n),indx(2,n), &
                      indx(3,n),indx(4,n), &
                      indx(5,n),indx(6,n) )
    end do

!-- remove far field on both scatter and main region
    call remove_farboth
!-- calculation on scatter region
    do n=1,SEIS_GEO*2
    call Lx[$1]_Ly[$2]_Lz[$3]( scaindx(1,n),scaindx(2,n), &
                      scaindx(3,n),scaindx(4,n), &
                      scaindx(5,n),scaindx(6,n) )
    end do
!-- recover true field in main domain
    call add_far2main

#ifdef CondFreeVHOC
!-- HOC velocity free surface
    if (freenode)  call Lx[$1]_Ly[$2]_Lz[$3]_VHOC
#endif

!-------------------------------------------------------------------------------
end subroutine in_Lx[$1]_Ly[$2]_Lz[$3]
!-------------------------------------------------------------------------------
])dnl
divert([0])dnl

_WRAPLXYZ(F,F,F)
_WRAPLXYZ(B,B,B)

_WRAPLXYZ(F,F,B)
_WRAPLXYZ(B,B,F)

_WRAPLXYZ(B,F,F)
_WRAPLXYZ(F,B,B)

_WRAPLXYZ(F,B,F)
_WRAPLXYZ(B,F,B)



!===============================================================================
!--- 3D bias FD of the MacCormack scheme in interior domain
!===============================================================================

divert([-1])
define([_LXYZ],
[dnl
!-------------------------------------------------------------------------------
subroutine Lx[$1]_Ly[$2]_Lz[$3](I1,I2,J1,J2,K1,K2)
!-------------------------------------------------------------------------------
    integer,intent(in) :: I1,I2,J1,J2,K1,K2
    integer :: i,j,k
    real(SP) :: DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz
    real(SP) :: DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz
    real(SP) :: DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz
    real(SP) :: rrho
    real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
    real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
    real(SP) :: tc45,tc46,tc56
#endif
    real(SP) :: E11,E22,E33,E12,E13,E23
! curvilinear
!$OMP PARALLEL DO DEFAULT(shared)  &
!$OMP PRIVATE(i,j,k, &
!$OMP   DxTxx,DxTxy,DxTxz,DxVx,DxVy,DxVz, &
!$OMP   DyTyy,DyTxy,DyTyz,DyVx,DyVy,DyVz, &
!$OMP   DzTzz,DzTxz,DzTyz,DzVx,DzVy,DzVz, &
!$OMP   lam,miu,lam2mu,rrho, &
!$OMP   E11,E22,E33,E12,E13,E23 )
    do k=K1,K2
    do j=J1,J2
    do i=I1,I2
       DxTxx = (              &
         m3d_FDx[$1]1(Txx,i,j,k) &
         m3d_FDx[$1]2(Txx,i,j,k) &
         )*xi_x(i)
       DxTxy = (              &
         m3d_FDx[$1]1(Txy,i,j,k) &
         m3d_FDx[$1]2(Txy,i,j,k) &
         )*xi_x(i)
       DxTxz = (              &
         m3d_FDx[$1]1(Txz,i,j,k) &
         m3d_FDx[$1]2(Txz,i,j,k) &
         )*xi_x(i)
       DxVx =  (              &
         m3d_FDx[$1]1(Vx,i,j,k)  &
         m3d_FDx[$1]2(Vx,i,j,k)  &
         )*xi_x(i)
       DxVy = (               &
         m3d_FDx[$1]1(Vy,i,j,k)  &
         m3d_FDx[$1]2(Vy,i,j,k)  &
         )*xi_x(i)
       DxVz = (               &
         m3d_FDx[$1]1(Vz,i,j,k)  &
         m3d_FDx[$1]2(Vz,i,j,k)  &
         )*xi_x(i)
    
       DyTyy = (              &
         m3d_FDy[$2]1(Tyy,i,j,k) &
         m3d_FDy[$2]2(Tyy,i,j,k) &
         )*eta_y(j)
       DyTxy = (              &
         m3d_FDy[$2]1(Txy,i,j,k) &
         m3d_FDy[$2]2(Txy,i,j,k) &
         )*eta_y(j)
       DyTyz = (              &
         m3d_FDy[$2]1(Tyz,i,j,k) &
         m3d_FDy[$2]2(Tyz,i,j,k) &
         )*eta_y(j)
       DyVx = (               &
         m3d_FDy[$2]1(Vx,i,j,k)  &
         m3d_FDy[$2]2(Vx,i,j,k)  &
         )*eta_y(j)
       DyVy = (               &
         m3d_FDy[$2]1(Vy,i,j,k)  &
         m3d_FDy[$2]2(Vy,i,j,k)  &
         )*eta_y(j)
       DyVz = (               &
         m3d_FDy[$2]1(Vz,i,j,k)  &
         m3d_FDy[$2]2(Vz,i,j,k)  &
         )*eta_y(j)
    
       DzTzz = (              &
         m3d_FDz[$3]1(Tzz,i,j,k) &
         m3d_FDz[$3]2(Tzz,i,j,k) &
         )*zeta_z(k)
       DzTxz = (              &
         m3d_FDz[$3]1(Txz,i,j,k) &
         m3d_FDz[$3]2(Txz,i,j,k) &
         )*zeta_z(k)
       DzTyz = (              &
         m3d_FDz[$3]1(Tyz,i,j,k) &
         m3d_FDz[$3]2(Tyz,i,j,k) &
         )*zeta_z(k)
#ifdef CondFreeVLOW
       if (freenode .and. k==nk2-1) then
       DzVx = (               &
         m22_FDz[$3]1(Vx,i,j,k)  &
         m22_FDz[$3]2(Vx,i,j,k)  &
         )*zeta_z(k)
       DzVy = (               &
         m22_FDz[$3]1(Vy,i,j,k)  &
         m22_FDz[$3]2(Vy,i,j,k)  &
         )*zeta_z(k)
       DzVz = (               &
         m22_FDz[$3]1(Vz,i,j,k)  &
         m22_FDz[$3]2(Vz,i,j,k)  &
         )*zeta_z(k)
       elseif (freenode .and. k==nk2-2) then
       DzVx = (               &
         m24_FDz[$3]1(Vx,i,j,k)  &
         m24_FDz[$3]2(Vx,i,j,k)  &
         )*zeta_z(k)
       DzVy = (               &
         m24_FDz[$3]1(Vy,i,j,k)  &
         m24_FDz[$3]2(Vy,i,j,k)  &
         )*zeta_z(k)
       DzVz = (               &
         m24_FDz[$3]1(Vz,i,j,k)  &
         m24_FDz[$3]2(Vz,i,j,k)  &
         )*zeta_z(k)
       else
#endif
       DzVx = (               &
         m3d_FDz[$3]1(Vx,i,j,k)  &
         m3d_FDz[$3]2(Vx,i,j,k)  &
         )*zeta_z(k)
       DzVy = (               &
         m3d_FDz[$3]1(Vy,i,j,k)  &
         m3d_FDz[$3]2(Vy,i,j,k)  &
         )*zeta_z(k)
       DzVz = (               &
         m3d_FDz[$3]1(Vz,i,j,k)  &
         m3d_FDz[$3]2(Vz,i,j,k)  &
         )*zeta_z(k)
#ifdef CondFreeVLOW
       end if
#endif
    
#ifdef WATER
       if (C66(i,j,k)<=SEIS_ZERO) then
          DxTxy=0.0
          DxTxz=0.0
          DyTxy=0.0
          DyTyz=0.0
          DzTxz=0.0
          DzTyz=0.0
       end if
#endif
    
       rrho=1.0/rho(i,j,k)
    
       hVx(i,j,k)= rrho*( DxTxx/z(k)+DyTxy/z(k)/xsin(i)+DzTxz                    &
            +(3.0_SP*Txz(i,j,k)+Txx(i,j,k)*xcot(i)-Tyy(i,j,k)*xcot(i))/z(k) )
       hVy(i,j,k)= rrho*( DxTxy/z(k)+DyTyy/z(k)/xsin(i)+DzTyz                    &
            +(2.0_SP*Txy(i,j,k)*xcot(i)+3.0_SP*Tyz(i,j,k))/z(k) )
       hVz(i,j,k)= rrho*( DxTxz/z(k)+DyTyz/z(k)/xsin(i)+DzTzz                    &
            +(2.0_SP*Tzz(i,j,k)-Txx(i,j,k)-Tyy(i,j,k)+Txz(i,j,k)*xcot(i))/z(k) )
    
#if defined AnisGene
       tc11=C11(i,j,k)
       tc12=C12(i,j,k)
       tc13=C13(i,j,k)
       tc14=C14(i,j,k)
       tc15=C15(i,j,k)
       tc16=C16(i,j,k)
       tc22=C22(i,j,k)
       tc23=C23(i,j,k)
       tc24=C24(i,j,k)
       tc25=C25(i,j,k)
       tc26=C26(i,j,k)
       tc33=C33(i,j,k)
       tc34=C34(i,j,k)
       tc35=C35(i,j,k)
       tc36=C36(i,j,k)
       tc44=C44(i,j,k)
       tc45=C45(i,j,k)
       tc46=C46(i,j,k)
       tc55=C55(i,j,k)
       tc56=C56(i,j,k)
       tc66=C66(i,j,k)
#elif defined AnisVTI
       tc11=C11(i,j,k)
       tc13=C13(i,j,k)
       tc33=C33(i,j,k)
       tc44=C44(i,j,k)
       tc66=C66(i,j,k)
    
       tc12=tc11-2.0_SP*tc66
       tc22=tc11
       tc23=tc13
       tc55=tc44
#else
       tc13=C13(i,j,k)
       tc66=C66(i,j,k)
    
       tc11=tc13+2.0_SP*tc66
       tc12=tc13
       tc22=tc11
       tc23=tc13
       tc33=tc11
       tc44=tc66
       tc55=tc66
#endif
    
       E11=(DxVx+Vz(i,j,k))/z(k)
       E22=(Vx(i,j,k)*xcot(i)+DyVy/xsin(i)+Vz(i,j,k))/z(k)
       E12=(DxVy+DyVx/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)
    
       E33=DzVz
       E13=(DxVz/z(k)+DzVx-Vx(i,j,k)/z(k))
       E23=(DyVz/z(k)/xsin(i)+DzVy-Vy(i,j,k)/z(k))
    
#if ! ( defined CondFreeCharac || defined CondFreeVHOC )
       if (freenode .and. k==nk2) then
#ifdef AnisGene
          E33= matEh2Ev(1,1,i,j)*E11+matEh2Ev(1,2,i,j)*E22+matEh2Ev(1,3,i,j)*E12 &
              +matF2Ev(1,1,i,j)*VzSrc(i,j), &
              +matF2Ev(1,2,i,j)*VySrc(i,j), &
              +matF2Ev(1,3,i,j)*VxSrc(i,j)
          E23= matEh2Ev(2,1,i,j)*E11+matEh2Ev(2,2,i,j)*E22+matEh2Ev(2,3,i,j)*E12 &
              +matF2Ev(2,1,i,j)*VzSrc(i,j), &
              +matF2Ev(2,2,i,j)*VySrc(i,j), &
              +matF2Ev(2,3,i,j)*VxSrc(i,j)
          E13= matEh2Ev(3,1,i,j)*E11+matEh2Ev(3,2,i,j)*E22+matEh2Ev(3,3,i,j)*E12 &
              +matF2Ev(3,1,i,j)*VzSrc(i,j), &
              +matF2Ev(3,2,i,j)*VySrc(i,j), &
              +matF2Ev(3,3,i,j)*VxSrc(i,j)
#else
          E33=(-tc13*E11-tc23*E22+VzSrc(i,j))/tc33
          if (tc44 > SEIS_ZERO) then !-- avoid devide zero in water
             E23=VySrc(i,j)/tc44
          else
             E23=0.0_SP
          end if
          if (tc55 > SEIS_ZERO) then
             E13=VxSrc(i,j)/tc55
          else
             E13=0.0_SP
          end if
#endif
       end if
#endif
    
#ifdef AnisGene
       hTxx(i,j,k)=tc11*E11+tc12*E22+tc13*E33+tc14*E23+tc15*E13+tc16*E12
       hTyy(i,j,k)=tc12*E11+tc22*E22+tc23*E33+tc24*E23+tc25*E13+tc26*E12
       hTzz(i,j,k)=tc13*E11+tc23*E22+tc33*E33+tc34*E23+tc35*E13+tc36*E12
       hTyz(i,j,k)=tc14*E11+tc24*E22+tc34*E33+tc44*E23+tc45*E13+tc46*E12
       hTxz(i,j,k)=tc15*E11+tc25*E22+tc35*E33+tc45*E23+tc55*E13+tc56*E12
       hTxy(i,j,k)=tc16*E11+tc26*E22+tc36*E33+tc46*E23+tc56*E13+tc66*E12
#else
       hTxx(i,j,k)=tc11*E11+tc12*E22+tc13*E33
       hTyy(i,j,k)=tc12*E11+tc22*E22+tc23*E33
       hTzz(i,j,k)=tc13*E11+tc23*E22+tc33*E33
       hTyz(i,j,k)=tc44*E23
       hTxz(i,j,k)=         tc55*E13
       hTxy(i,j,k)=                  tc66*E12
#endif
    
    end do
    end do
    end do
!$OMP END PARALLEL DO
!-------------------------------------------------------------------------------
end subroutine Lx[$1]_Ly[$2]_Lz[$3]
!-------------------------------------------------------------------------------
])dnl
divert([0])dnl

_LXYZ(F,F,F)
_LXYZ(B,B,B)

_LXYZ(F,F,B)
_LXYZ(B,B,F)

_LXYZ(B,F,F)
_LXYZ(F,B,B)

_LXYZ(F,B,F)
_LXYZ(B,F,B)


!===============================================================================
!--  use Compact MacCormack scheme to calculate velocities fd with
!--  respect to eta and assemble the right hand side to update stresses
!===============================================================================

#ifdef CondFreeVHOC

divert([-1])
define([_LXYZVHOC_ZF],
[dnl
!-------------------------------------------------------------------------------
subroutine Lx[$1]_Ly[$2]_LzF_VHOC
!-------------------------------------------------------------------------------
    integer :: i,j,k,n
    real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
    real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
    real(SP) :: tc45,tc46,tc56
#endif
    real(SP) :: rhs_Dz,lhs_Dz
    real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                                 DyVx,DyVy,DyVz, &
                                 DzVx,DzVy,DzVz
    real(SP) :: E11,E22,E33,E12,E13,E23
    
    if (.not. freenode) return
    
!-- velocity fd --
    loop_eta: do j=nj1,nj2
    loop_xi:  do i=ni1,ni2
       !-- k=nz --
       n=LenFD+1; k=nk2
    
#if ! defined AnisGene
#if defined AnisVTI
       tc13=C13(i,j,k)
       tc33=C33(i,j,k)
       tc44=C44(i,j,k)
    
       tc23=tc13
       tc55=tc44
#else
       tc13=C13(i,j,k)
       tc66=C66(i,j,k)
    
       tc11=tc13+2.0_SP*tc66
       tc23=tc13
       tc33=tc11
       tc44=tc66
       tc55=tc66
#endif
#endif
    
       DxVx(n) =  (          &
         m3d_FDx[$1]1(Vx,i,j,k) &
         m3d_FDx[$1]2(Vx,i,j,k) &
         )*xi_x(i)
       DxVy(n) = (           &
         m3d_FDx[$1]1(Vy,i,j,k) &
         m3d_FDx[$1]2(Vy,i,j,k) &
         )*xi_x(i)
       DxVz(n) = (           &
         m3d_FDx[$1]1(Vz,i,j,k) &
         m3d_FDx[$1]2(Vz,i,j,k) &
         )*xi_x(i)
    
       DyVx(n) =  (          &
         m3d_FDy[$2]1(Vx,i,j,k) &
         m3d_FDy[$2]2(Vx,i,j,k) &
         )*eta_y(j)
       DyVy(n) = (           &
         m3d_FDy[$2]1(Vy,i,j,k) &
         m3d_FDy[$2]2(Vy,i,j,k) &
         )*eta_y(j)
       DyVz(n) = (           &
         m3d_FDy[$2]1(Vz,i,j,k) &
         m3d_FDy[$2]2(Vz,i,j,k) &
         )*eta_y(j)
    
       E11=(DxVx(n)+Vz(i,j,k))/z(k)
       E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
    
#ifdef AnisGene
       E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)
       E33= matEh2Ev(1,1,i,j)*E11+matEh2Ev(1,2,i,j)*E22+matEh2Ev(1,3,i,j)*E12 &
              +matF2Ev(1,1,i,j)*VzSrc(i,j), &
              +matF2Ev(1,2,i,j)*VySrc(i,j), &
              +matF2Ev(1,3,i,j)*VxSrc(i,j)
       E23= matEh2Ev(2,1,i,j)*E11+matEh2Ev(2,2,i,j)*E22+matEh2Ev(2,3,i,j)*E12 &
              +matF2Ev(2,1,i,j)*VzSrc(i,j), &
              +matF2Ev(2,2,i,j)*VySrc(i,j), &
              +matF2Ev(2,3,i,j)*VxSrc(i,j)
       E13= matEh2Ev(3,1,i,j)*E11+matEh2Ev(3,2,i,j)*E22+matEh2Ev(3,3,i,j)*E12 &
              +matF2Ev(3,1,i,j)*VzSrc(i,j), &
              +matF2Ev(3,2,i,j)*VySrc(i,j), &
              +matF2Ev(3,3,i,j)*VxSrc(i,j)
#else
       E33=(-tc13*E11-tc23*E22+VzSrc(i,j))/tc33
       E23=VySrc(i,j)/tc44
       E13=VxSrc(i,j)/tc55
#endif
    
       DzVz(n) = E33/zeta_z(k)
       DzVy(n) = (E23-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k))/zeta_z(k)
       DzVx(n) = (E13-DxVz(n)/z(k)+Vx(i,j,k)/z(k))/zeta_z(k)
    
       !-- j=ny-LenFD+1:ny-1 --
    do n=LenFD,2,-1
       k=nk2-LenFD +n-1
    
       DxVx(n) =  (          &
         m3d_FDx[$1]1(Vx,i,j,k) &
         m3d_FDx[$1]2(Vx,i,j,k) &
         )*xi_x(i)
       DxVy(n) = (           &
         m3d_FDx[$1]1(Vy,i,j,k) &
         m3d_FDx[$1]2(Vy,i,j,k) &
         )*xi_x(i)
       DxVz(n) = (           &
         m3d_FDx[$1]1(Vz,i,j,k) &
         m3d_FDx[$1]2(Vz,i,j,k) &
         )*xi_x(i)
    
       DyVx(n) =  (          &
         m3d_FDy[$2]1(Vx,i,j,k) &
         m3d_FDy[$2]2(Vx,i,j,k) &
         )*eta_y(j)
       DyVy(n) = (           &
         m3d_FDy[$2]1(Vy,i,j,k) &
         m3d_FDy[$2]2(Vy,i,j,k) &
         )*eta_y(j)
       DyVz(n) = (           &
         m3d_FDy[$2]1(Vz,i,j,k) &
         m3d_FDy[$2]2(Vz,i,j,k) &
         )*eta_y(j)
    
       rhs_Dz= (                  &
         m3d_HOCzF1_RHS(Vx,i,j,k) &
         m3d_HOCzF2_RHS(Vx,i,j,k) &
         )
       lhs_Dz= (                  &
         vec_HOC_F_LHS(DzVx,n)    &
         )
       DzVx(n)=rhs_Dz - lhs_Dz
    
       rhs_Dz= (                  &
         m3d_HOCzF1_RHS(Vy,i,j,k) &
         m3d_HOCzF2_RHS(Vy,i,j,k) &
         )
       lhs_Dz= (                  &
         vec_HOC_F_LHS(DzVy,n)    &
         )
       DzVy(n)=rhs_Dz - lhs_Dz
    
       rhs_Dz= (                  &
         m3d_HOCzF1_RHS(Vz,i,j,k) &
         m3d_HOCzF2_RHS(Vz,i,j,k) &
         )
       lhs_Dz= (                  &
         vec_HOC_F_LHS(DzVz,n)    &
         )
       DzVz(n)=rhs_Dz - lhs_Dz
    
    end do
    do n=2,LenFD+1
       k=nk2-LenFD +n-1
       DzVx(n)=DzVx(n)*zeta_z(k)
       DzVy(n)=DzVy(n)*zeta_z(k)
       DzVz(n)=DzVz(n)*zeta_z(k)
    
#if defined AnisGene
       tc11=C11(i,j,k)
       tc12=C12(i,j,k)
       tc13=C13(i,j,k)
       tc14=C14(i,j,k)
       tc15=C15(i,j,k)
       tc16=C16(i,j,k)
       tc22=C22(i,j,k)
       tc23=C23(i,j,k)
       tc24=C24(i,j,k)
       tc25=C25(i,j,k)
       tc26=C26(i,j,k)
       tc33=C33(i,j,k)
       tc34=C34(i,j,k)
       tc35=C35(i,j,k)
       tc36=C36(i,j,k)
       tc44=C44(i,j,k)
       tc45=C45(i,j,k)
       tc46=C46(i,j,k)
       tc55=C55(i,j,k)
       tc56=C56(i,j,k)
       tc66=C66(i,j,k)
#elif defined AnisVTI
       tc11=C11(i,j,k)
       tc13=C13(i,j,k)
       tc33=C33(i,j,k)
       tc44=C44(i,j,k)
       tc66=C66(i,j,k)
    
       tc12=tc11-2.0_SP*tc66
       tc22=tc11
       tc23=tc13
       tc55=tc44
#else
       tc13=C13(i,j,k)
       tc66=C66(i,j,k)
    
       tc11=tc13+2.0_SP*tc66
       tc12=tc13
       tc22=tc11
       tc23=tc13
       tc33=tc11
       tc44=tc66
       tc55=tc66
#endif
    
       E11=(DxVx(n)+Vz(i,j,k))/z(k)
       E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
       E33=DzVz(n)
       E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)
       E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))
       E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))
    
#ifdef AnisGene
       hTxx(i,j,k)=tc11*E11+tc12*E22+tc13*E33+tc14*E23+tc15*E13+tc16*E12
       hTyy(i,j,k)=tc12*E11+tc22*E22+tc23*E33+tc24*E23+tc25*E13+tc26*E12
       hTzz(i,j,k)=tc13*E11+tc23*E22+tc33*E33+tc34*E23+tc35*E13+tc36*E12
       hTyz(i,j,k)=tc14*E11+tc24*E22+tc34*E33+tc44*E23+tc45*E13+tc46*E12
       hTxz(i,j,k)=tc15*E11+tc25*E22+tc35*E33+tc45*E23+tc55*E13+tc56*E12
       hTxy(i,j,k)=tc16*E11+tc26*E22+tc36*E33+tc46*E23+tc56*E13+tc66*E12
#else
       hTxx(i,j,k)=tc11*E11+tc12*E22+tc13*E33
       hTyy(i,j,k)=tc12*E11+tc22*E22+tc23*E33
       hTzz(i,j,k)=tc13*E11+tc23*E22+tc33*E33
       hTyz(i,j,k)=tc44*E23
       hTxz(i,j,k)=         tc55*E13
       hTxy(i,j,k)=                  tc66*E12
#endif
       if (k==nk2) then
          hTxz(i,j,k)=VxSrc(i,j)
          hTyz(i,j,k)=VySrc(i,j)
          hTzz(i,j,k)=VzSrc(i,j)
       end if
    end do
    end do loop_xi
    end do loop_eta
!-------------------------------------------------------------------------------
end subroutine Lx[$1]_Ly[$2]_LzF_VHOC
!-------------------------------------------------------------------------------
])dnl
divert([0])dnl


divert([-1])
define([_LXYZVHOC_ZB],
[dnl
!-------------------------------------------------------------------------------
subroutine Lx[$1]_Ly[$2]_LzB_VHOC
!-------------------------------------------------------------------------------
    integer :: i,j,k,n
    real(SP) :: tc11,tc12,tc13,tc22,tc23,tc33,tc44,tc55,tc66
#ifdef AnisGene
    real(SP) :: tc14,tc15,tc16,tc24,tc25,tc26,tc34,tc35,tc36
    real(SP) :: tc45,tc46,tc56
#endif
    real(SP) :: rhs_Dz,lhs_Dz
    real(SP),dimension(1:LenFD+1) :: DxVx,DxVy,DxVz, &
                                 DyVx,DyVy,DyVz, &
                                 DzVx,DzVy,DzVz
    real(SP) :: E11,E22,E33,E12,E13,E23

    if (.not. freenode) return

!-- velocity fd --
    loop_eta: do j=nj1,nj2
    loop_xi:  do i=ni1,ni2
       !-- k=nz-LenFD --
       n=1; k=nk2-LenFD
       DzVx(n) =  (          &
         m3d_FDzB1(Vx,i,j,k) &
         m3d_FDzB2(Vx,i,j,k) &
         )
       DzVy(n) = (           &
         m3d_FDzB1(Vy,i,j,k) &
         m3d_FDzB2(Vy,i,j,k) &
         )
       DzVz(n) = (           &
         m3d_FDzB1(Vz,i,j,k) &
         m3d_FDzB2(Vz,i,j,k) &
         )
    
       !-- k=nk2 --
       n=LenFD+1; k=nk2
       DxVx(n) =  (          &
         m3d_FDx[$1]1(Vx,i,j,k) &
         m3d_FDx[$1]2(Vx,i,j,k) &
         )*xi_x(i)
       DxVy(n) = (           &
         m3d_FDx[$1]1(Vy,i,j,k) &
         m3d_FDx[$1]2(Vy,i,j,k) &
         )*xi_x(i)
       DxVz(n) = (           &
         m3d_FDx[$1]1(Vz,i,j,k) &
         m3d_FDx[$1]2(Vz,i,j,k) &
         )*xi_x(i)
    
       DyVx(n) =  (          &
         m3d_FDy[$2]1(Vx,i,j,k) &
         m3d_FDy[$2]2(Vx,i,j,k) &
         )*eta_y(j)
       DyVy(n) = (           &
         m3d_FDy[$2]1(Vy,i,j,k) &
         m3d_FDy[$2]2(Vy,i,j,k) &
         )*eta_y(j)
       DyVz(n) = (           &
         m3d_FDy[$2]1(Vz,i,j,k) &
         m3d_FDy[$2]2(Vz,i,j,k) &
         )*eta_y(j)
    
#if ! defined AnisGene
#if defined AnisVTI
       tc13=C13(i,j,k)
       tc33=C33(i,j,k)
       tc44=C44(i,j,k)
    
       tc23=tc13
       tc55=tc44
#else
       tc13=C13(i,j,k)
       tc66=C66(i,j,k)
    
       tc11=tc13+2.0_SP*tc66
       tc23=tc13
       tc33=tc11
       tc44=tc66
       tc55=tc66
#endif
#endif
    
       E11=(DxVx(n)+Vz(i,j,k))/z(k)
       E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
    
#ifdef AnisGene
       E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)
       E33= matEh2Ev(1,1,i,j)*E11+matEh2Ev(1,2,i,j)*E22+matEh2Ev(1,3,i,j)*E12 &
              +matF2Ev(1,1,i,j)*VzSrc(i,j), &
              +matF2Ev(1,2,i,j)*VySrc(i,j), &
              +matF2Ev(1,3,i,j)*VxSrc(i,j)
       E23= matEh2Ev(2,1,i,j)*E11+matEh2Ev(2,2,i,j)*E22+matEh2Ev(2,3,i,j)*E12 &
              +matF2Ev(2,1,i,j)*VzSrc(i,j), &
              +matF2Ev(2,2,i,j)*VySrc(i,j), &
              +matF2Ev(2,3,i,j)*VxSrc(i,j)
       E13= matEh2Ev(3,1,i,j)*E11+matEh2Ev(3,2,i,j)*E22+matEh2Ev(3,3,i,j)*E12 &
              +matF2Ev(3,1,i,j)*VzSrc(i,j), &
              +matF2Ev(3,2,i,j)*VySrc(i,j), &
              +matF2Ev(3,3,i,j)*VxSrc(i,j)
#else
       E33=(-tc13*E11-tc23*E22+VzSrc(i,j))/tc33
       E23=VySrc(i,j)/tc44
       E13=VxSrc(i,j)/tc55
#endif
    
       DzVz(n) = E33/zeta_z(k)
       DzVy(n) = (E23-DyVz(n)/z(k)/xsin(i)+Vy(i,j,k)/z(k))/zeta_z(k)
       DzVx(n) = (E13-DxVz(n)/z(k)+Vx(i,j,k)/z(k))/zeta_z(k)
    
       !-- k=nk2-LenFD+1:nk2-1 --
    do n=2,LenFD
       k=nk2-LenFD +n-1
    
       DxVx(n) =  (          &
         m3d_FDx[$1]1(Vx,i,j,k) &
         m3d_FDx[$1]2(Vx,i,j,k) &
         )*xi_x(i)
       DxVy(n) = (           &
         m3d_FDx[$1]1(Vy,i,j,k) &
         m3d_FDx[$1]2(Vy,i,j,k) &
         )*xi_x(i)
       DxVz(n) = (           &
         m3d_FDx[$1]1(Vz,i,j,k) &
         m3d_FDx[$1]2(Vz,i,j,k) &
         )*xi_x(i)
    
       DyVx(n) =  (          &
         m3d_FDy[$2]1(Vx,i,j,k) &
         m3d_FDy[$2]2(Vx,i,j,k) &
         )*eta_y(j)
       DyVy(n) = (           &
         m3d_FDy[$2]1(Vy,i,j,k) &
         m3d_FDy[$2]2(Vy,i,j,k) &
         )*eta_y(j)
       DyVz(n) = (           &
         m3d_FDy[$2]1(Vz,i,j,k) &
         m3d_FDy[$2]2(Vz,i,j,k) &
         )*eta_y(j)
    
       rhs_Dz= (                  &
         m3d_HOCzB1_RHS(Vx,i,j,k) &
         m3d_HOCzB2_RHS(Vx,i,j,k) &
         )
       lhs_Dz= (                  &
         vec_HOC_B_LHS(DzVx,n)    &
         )
       DzVx(n)=rhs_Dz - lhs_Dz
    
       rhs_Dz= (                  &
         m3d_HOCzB1_RHS(Vy,i,j,k) &
         m3d_HOCzB2_RHS(Vy,i,j,k) &
         )
       lhs_Dz= (                  &
         vec_HOC_B_LHS(DzVy,n)    &
         )
       DzVy(n)=rhs_Dz - lhs_Dz
    
       rhs_Dz= (                  &
         m3d_HOCzB1_RHS(Vz,i,j,k) &
         m3d_HOCzB2_RHS(Vz,i,j,k) &
         )
       lhs_Dz= (                  &
         vec_HOC_B_LHS(DzVz,n)    &
         )
       DzVz(n)=rhs_Dz - lhs_Dz
    end do
    
    do n=2,LenFD+1
       k=nk2-LenFD +n-1
       DzVx(n)=DzVx(n)*zeta_z(k)
       DzVy(n)=DzVy(n)*zeta_z(k)
       DzVz(n)=DzVz(n)*zeta_z(k)
    
#if defined AnisGene
       tc11=C11(i,j,k)
       tc12=C12(i,j,k)
       tc13=C13(i,j,k)
       tc14=C14(i,j,k)
       tc15=C15(i,j,k)
       tc16=C16(i,j,k)
       tc22=C22(i,j,k)
       tc23=C23(i,j,k)
       tc24=C24(i,j,k)
       tc25=C25(i,j,k)
       tc26=C26(i,j,k)
       tc33=C33(i,j,k)
       tc34=C34(i,j,k)
       tc35=C35(i,j,k)
       tc36=C36(i,j,k)
       tc44=C44(i,j,k)
       tc45=C45(i,j,k)
       tc46=C46(i,j,k)
       tc55=C55(i,j,k)
       tc56=C56(i,j,k)
       tc66=C66(i,j,k)
#elif defined AnisVTI
       tc11=C11(i,j,k)
       tc13=C13(i,j,k)
       tc33=C33(i,j,k)
       tc44=C44(i,j,k)
       tc66=C66(i,j,k)
    
       tc12=tc11-2.0_SP*tc66
       tc22=tc11
       tc23=tc13
       tc55=tc44
#else
       tc13=C13(i,j,k)
       tc66=C66(i,j,k)
    
       tc11=tc13+2.0_SP*tc66
       tc12=tc13
       tc22=tc11
       tc23=tc13
       tc33=tc11
       tc44=tc66
       tc55=tc66
#endif
    
       E11=(DxVx(n)+Vz(i,j,k))/z(k)
       E22=(Vx(i,j,k)*xcot(i)+DyVy(n)/xsin(i)+Vz(i,j,k))/z(k)
       E33=DzVz(n)
       E12=(DxVy(n)+DyVx(n)/xsin(i)-Vy(i,j,k)*xcot(i))/z(k)
       E13=(DxVz(n)/z(k)+DzVx(n)-Vx(i,j,k)/z(k))
       E23=(DyVz(n)/z(k)/xsin(i)+DzVy(n)-Vy(i,j,k)/z(k))
    
#ifdef AnisGene
       hTxx(i,j,k)=tc11*E11+tc12*E22+tc13*E33+tc14*E23+tc15*E13+tc16*E12
       hTyy(i,j,k)=tc12*E11+tc22*E22+tc23*E33+tc24*E23+tc25*E13+tc26*E12
       hTzz(i,j,k)=tc13*E11+tc23*E22+tc33*E33+tc34*E23+tc35*E13+tc36*E12
       hTyz(i,j,k)=tc14*E11+tc24*E22+tc34*E33+tc44*E23+tc45*E13+tc46*E12
       hTxz(i,j,k)=tc15*E11+tc25*E22+tc35*E33+tc45*E23+tc55*E13+tc56*E12
       hTxy(i,j,k)=tc16*E11+tc26*E22+tc36*E33+tc46*E23+tc56*E13+tc66*E12
#else
       hTxx(i,j,k)=tc11*E11+tc12*E22+tc13*E33
       hTyy(i,j,k)=tc12*E11+tc22*E22+tc23*E33
       hTzz(i,j,k)=tc13*E11+tc23*E22+tc33*E33
       hTyz(i,j,k)=tc44*E23
       hTxz(i,j,k)=         tc55*E13
       hTxy(i,j,k)=                  tc66*E12
#endif
       if (k==nk2) then
          hTxz(i,j,k)=VxSrc(i,j)
          hTyz(i,j,k)=VySrc(i,j)
          hTzz(i,j,k)=VzSrc(i,j)
       end if
    end do
    end do loop_xi
    end do loop_eta
!-------------------------------------------------------------------------------
end subroutine Lx[$1]_Ly[$2]_LzB_VHOC
!-------------------------------------------------------------------------------
])dnl
divert([0])dnl


_LXYZVHOC_ZF(F,F)
_LXYZVHOC_ZB(B,B)

_LXYZVHOC_ZB(F,F)
_LXYZVHOC_ZF(B,B)

_LXYZVHOC_ZF(B,F)
_LXYZVHOC_ZB(F,B)

_LXYZVHOC_ZF(F,B)
_LXYZVHOC_ZB(B,F)

#endif

!===============================================================================
!--- Far-Near field coupling correction
!===============================================================================

!-------------------------------------------------------------------------------
subroutine read_farfield
!-------------------------------------------------------------------------------
!-- check if coupling is needed
    if ( .not. flag_couple_node ) return
    if ( cur_nt > ntwin_Far_field ) return

!-- read Far data
    read(fid_Far) Vx_Far_1D,Vy_Far_1D , Vz_Far_1D ,   &
                  Txx_Far_1D, Tyy_Far_1D, Tzz_Far_1D, &
                  Tyz_Far_1D, Txz_Far_1D, Txy_Far_1D 

!-------------------------------------------------------------------------------
end subroutine read_farfield
!-------------------------------------------------------------------------------

!!-------------------------------------------------------------------------------
!subroutine add_far2scat_main
!!-------------------------------------------------------------------------------
!    integer :: i,j,k,m,n
!
!!-- check if coupling is needed
!    if ( .not. flag_couple_node ) return
!    if ( cur_nt > ntwin_Far_field ) return
!
!    m=0
!    do n=1,2*SEIS_GEO
!       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
!       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
!       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
!          m=m+1
!          if (     (n==1 .and. i<indx_couple_grid_center(1,1)              &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==2 .and. i>indx_couple_grid_center(2,1)              &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==3 .and. j<indx_couple_grid_center(1,2)              &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==4 .and. j>indx_couple_grid_center(2,2)              &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==5 .and. k<indx_couple_grid_center(1,3)              &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7) ) &
!              .or. (n==6 .and. k>indx_couple_grid_center(2,3)              &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7) ) &
!               ) then
!               Vx(i,j,k)= Vx(i,j,k)+ Vx_Far_1D(m)
!               Vy(i,j,k)= Vy(i,j,k)+ Vy_Far_1D(m)
!               Vz(i,j,k)= Vz(i,j,k)+ Vz_Far_1D(m)
!              Txx(i,j,k)=Txx(i,j,k)+Txx_Far_1D(m)
!              Tyy(i,j,k)=Tyy(i,j,k)+Tyy_Far_1D(m)
!              Tzz(i,j,k)=Tzz(i,j,k)+Tzz_Far_1D(m)
!              Tyz(i,j,k)=Tyz(i,j,k)+Tyz_Far_1D(m)
!              Txz(i,j,k)=Txz(i,j,k)+Txz_Far_1D(m)
!              Txy(i,j,k)=Txy(i,j,k)+Txy_Far_1D(m)
!          end if
!       end do
!       end do
!       end do
!    end do
!!-------------------------------------------------------------------------------
!end subroutine add_far2scat_main
!!-------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------------
!subroutine add_far2scat_boundary
!!-------------------------------------------------------------------------------
!    integer :: i,j,k,m,n
!
!!-- check if coupling is needed
!    if ( .not. flag_couple_node ) return
!    if ( cur_nt > ntwin_Far_field ) return
!
!    m=0
!    do n=1,2*SEIS_GEO
!       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
!       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
!       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
!          m=m+1
!          if (     (n==1 .and. i<indx_couple_grid_center(1,1)             &
!                         .and. (indx(3,7)>j .or. j>indx(4,7)    &
!                           .or. indx(5,7)>k .or. k>indx(6,7)) ) &
!              .or. (n==2 .and. i>indx_couple_grid_center(2,1)             &
!                         .and. (indx(3,7)>j .or. j>indx(4,7)    &
!                          .or.  indx(5,7)>k .or. k>indx(6,7)) ) &
!              .or. (n==3 .and. j<indx_couple_grid_center(1,2)             &
!                         .and. (indx(1,7)>i .or. i>indx(2,7)    &
!                          .or.  indx(5,7)>k .or. k>indx(6,7)) ) &
!              .or. (n==4 .and. j>indx_couple_grid_center(2,2)             &
!                         .and. (indx(1,7)>i .or. i>indx(2,7)    &
!                          .or.  indx(5,7)>k .or. k>indx(6,7)) ) &
!              .or. (n==5 .and. k<indx_couple_grid_center(1,3)             &
!                         .and. (indx(3,7)>j .or. j>indx(4,7)    &
!                          .or.  indx(1,7)>i .or. i>indx(2,7)) ) &
!              .or. (n==6 .and. k>indx_couple_grid_center(2,3)             &
!                         .and. (indx(3,7)>j .or. j>indx(4,7)    &
!                          .or.  indx(1,7)>i .or. i>indx(2,7)) ) &
!               ) then
!               Vx(i,j,k)= Vx(i,j,k)+ Vx_Far_1D(m)
!               Vy(i,j,k)= Vy(i,j,k)+ Vy_Far_1D(m)
!               Vz(i,j,k)= Vz(i,j,k)+ Vz_Far_1D(m)
!              Txx(i,j,k)=Txx(i,j,k)+Txx_Far_1D(m)
!              Tyy(i,j,k)=Tyy(i,j,k)+Tyy_Far_1D(m)
!              Tzz(i,j,k)=Tzz(i,j,k)+Tzz_Far_1D(m)
!              Tyz(i,j,k)=Tyz(i,j,k)+Tyz_Far_1D(m)
!              Txz(i,j,k)=Txz(i,j,k)+Txz_Far_1D(m)
!              Txy(i,j,k)=Txy(i,j,k)+Txy_Far_1D(m)
!          end if
!       end do
!       end do
!       end do
!    end do
!!-------------------------------------------------------------------------------
!end subroutine add_far2scat_boundary
!!-------------------------------------------------------------------------------
!
!!-------------------------------------------------------------------------------
!subroutine add_far2main
!!-------------------------------------------------------------------------------
!    integer :: i,j,k,m,n
!
!!-- check if coupling is needed
!    if ( .not. flag_couple_node ) return
!    if ( cur_nt > ntwin_Far_field ) return
!
!    m=0
!    do n=1,2*SEIS_GEO
!       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
!       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
!       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
!          m=m+1
!          !if (     (n==1 .and. i>=indx_couple_grid_center(1,1)) &
!          !    .or. (n==2 .and. i<=indx_couple_grid_center(2,1)) &
!          !    .or. (n==3 .and. j>=indx_couple_grid_center(1,2)) &
!          !    .or. (n==4 .and. j<=indx_couple_grid_center(2,2)) &
!          !    .or. (n==5 .and. k>=indx_couple_grid_center(1,3)) &
!          !    .or. (n==6 .and. k<=indx_couple_grid_center(2,3)) ) then
!          if (     (n==1 .and. i>=indx_couple_grid_center(1,1)   &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==2 .and. i<=indx_couple_grid_center(2,1)   &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==3 .and. j>=indx_couple_grid_center(1,2)   &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==4 .and. j<=indx_couple_grid_center(2,2)   &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7)   &
!                         .and. indx(5,7)<=k .and. k<=indx(6,7) ) &
!              .or. (n==5 .and. k>=indx_couple_grid_center(1,3)   &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7) ) &
!              .or. (n==6 .and. k<=indx_couple_grid_center(2,3)   &
!                         .and. indx(3,7)<=j .and. j<=indx(4,7)   &
!                         .and. indx(1,7)<=i .and. i<=indx(2,7) ) &
!               ) then
!               Vx(i,j,k)= Vx(i,j,k)+ Vx_Far_1D(m)
!               Vy(i,j,k)= Vy(i,j,k)+ Vy_Far_1D(m)
!               Vz(i,j,k)= Vz(i,j,k)+ Vz_Far_1D(m)
!              Txx(i,j,k)=Txx(i,j,k)+Txx_Far_1D(m)
!              Tyy(i,j,k)=Tyy(i,j,k)+Tyy_Far_1D(m)
!              Tzz(i,j,k)=Tzz(i,j,k)+Tzz_Far_1D(m)
!              Tyz(i,j,k)=Tyz(i,j,k)+Tyz_Far_1D(m)
!              Txz(i,j,k)=Txz(i,j,k)+Txz_Far_1D(m)
!              Txy(i,j,k)=Txy(i,j,k)+Txy_Far_1D(m)
!          end if
!       end do
!       end do
!       end do
!    end do
!!-------------------------------------------------------------------------------
!end subroutine add_far2main
!!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine add_far2scat_main
!-------------------------------------------------------------------------------
    integer :: i,j,k,m,n

!-- check if coupling is needed
    if ( .not. flag_couple_node ) return
    if ( cur_nt > ntwin_Far_field ) return

    m=0
    do n=1,2*SEIS_GEO
       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
          m=m+1
          if (  (      i<indx_couple_grid_center(1,1)              &
                  .or. i>indx_couple_grid_center(2,1)              &
                  .or. j<indx_couple_grid_center(1,2)              &
                  .or. j>indx_couple_grid_center(2,2)              &
                  .or. k<indx_couple_grid_center(1,3)              &
                  .or. k>indx_couple_grid_center(2,3)              &
                 )  &  !- in scattering region
              .and.                &
                 (      i>=ni1+LenFD &
                  .and. i<=ni2-LenFD &
                  .and. j>=nj1+LenFD &
                  .and. j<=nj2-LenFD &
                  .and. k>=nk1+LenFD &
                  .and. k<=nk2-LenFD &
                  ) &  !- but not in MPI layer
               ) then
               Vx(i,j,k)= Vx(i,j,k)+ Vx_Far_1D(m)
               Vy(i,j,k)= Vy(i,j,k)+ Vy_Far_1D(m)
               Vz(i,j,k)= Vz(i,j,k)+ Vz_Far_1D(m)
              Txx(i,j,k)=Txx(i,j,k)+Txx_Far_1D(m)
              Tyy(i,j,k)=Tyy(i,j,k)+Tyy_Far_1D(m)
              Tzz(i,j,k)=Tzz(i,j,k)+Tzz_Far_1D(m)
              Tyz(i,j,k)=Tyz(i,j,k)+Tyz_Far_1D(m)
              Txz(i,j,k)=Txz(i,j,k)+Txz_Far_1D(m)
              Txy(i,j,k)=Txy(i,j,k)+Txy_Far_1D(m)
          end if
       end do
       end do
       end do
    end do
!-------------------------------------------------------------------------------
end subroutine add_far2scat_main
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine add_far2scat_boundary
!-------------------------------------------------------------------------------
    integer :: i,j,k,m,n

!-- check if coupling is needed
    if ( .not. flag_couple_node ) return
    if ( cur_nt > ntwin_Far_field ) return

    m=0
    do n=1,2*SEIS_GEO
       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
          m=m+1
          if (  (      i<indx_couple_grid_center(1,1)              &
                  .or. i>indx_couple_grid_center(2,1)              &
                  .or. j<indx_couple_grid_center(1,2)              &
                  .or. j>indx_couple_grid_center(2,2)              &
                  .or. k<indx_couple_grid_center(1,3)              &
                  .or. k>indx_couple_grid_center(2,3)              &
                 )  &  !- in scattering region
              .and.                &
                 (     i<ni1+LenFD &
                  .or. i>ni2-LenFD &
                  .or. j<nj1+LenFD &
                  .or. j>nj2-LenFD &
                  .or. k<nk1+LenFD &
                  .or. k>nk2-LenFD &
                  ) &  !- in MPI layer
               ) then
               Vx(i,j,k)= Vx(i,j,k)+ Vx_Far_1D(m)
               Vy(i,j,k)= Vy(i,j,k)+ Vy_Far_1D(m)
               Vz(i,j,k)= Vz(i,j,k)+ Vz_Far_1D(m)
              Txx(i,j,k)=Txx(i,j,k)+Txx_Far_1D(m)
              Tyy(i,j,k)=Tyy(i,j,k)+Tyy_Far_1D(m)
              Tzz(i,j,k)=Tzz(i,j,k)+Tzz_Far_1D(m)
              Tyz(i,j,k)=Tyz(i,j,k)+Tyz_Far_1D(m)
              Txz(i,j,k)=Txz(i,j,k)+Txz_Far_1D(m)
              Txy(i,j,k)=Txy(i,j,k)+Txy_Far_1D(m)
          end if
       end do
       end do
       end do
    end do
!-------------------------------------------------------------------------------
end subroutine add_far2scat_boundary
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine remove_farboth
!-------------------------------------------------------------------------------
    integer :: i,j,k,m,n

!-- check if coupling is needed
    if ( .not. flag_couple_node ) return
    if ( cur_nt > ntwin_Far_field ) return

    m=0
    do n=1,2*SEIS_GEO
       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
          m=m+1
          !if (     (n==1 .and. i<indx_couple_grid_center(1,1)) &
          !    .or. (n==2 .and. i>indx_couple_grid_center(2,1)) &
          !    .or. (n==3 .and. i<indx_couple_grid_center(1,2)) &
          !    .or. (n==4 .and. i>indx_couple_grid_center(2,2)) &
          !    .or. (n==5 .and. i<indx_couple_grid_center(1,3)) &
          !    .or. (n==6 .and. i>indx_couple_grid_center(2,3)) ) then
               Vx(i,j,k)= Vx(i,j,k)- Vx_Far_1D(m)
               Vy(i,j,k)= Vy(i,j,k)- Vy_Far_1D(m)
               Vz(i,j,k)= Vz(i,j,k)- Vz_Far_1D(m)
              Txx(i,j,k)=Txx(i,j,k)-Txx_Far_1D(m)
              Tyy(i,j,k)=Tyy(i,j,k)-Tyy_Far_1D(m)
              Tzz(i,j,k)=Tzz(i,j,k)-Tzz_Far_1D(m)
              Tyz(i,j,k)=Tyz(i,j,k)-Tyz_Far_1D(m)
              Txz(i,j,k)=Txz(i,j,k)-Txz_Far_1D(m)
              Txy(i,j,k)=Txy(i,j,k)-Txy_Far_1D(m)
          !end if
       end do
       end do
       end do
    end do
!-------------------------------------------------------------------------------
end subroutine remove_farboth
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine add_far2main
!-------------------------------------------------------------------------------
    integer :: i,j,k,m,n

!-- check if coupling is needed
    if ( .not. flag_couple_node ) return
    if ( cur_nt > ntwin_Far_field ) return

    m=0
    do n=1,2*SEIS_GEO
       do k=indx_Far_1D(5,n),indx_Far_1D(6,n)
       do j=indx_Far_1D(3,n),indx_Far_1D(4,n)
       do i=indx_Far_1D(1,n),indx_Far_1D(2,n)
          m=m+1
          if (      i>=indx_couple_grid_center(1,1)              &
              .and. i<=indx_couple_grid_center(2,1)              &
              .and. j>=indx_couple_grid_center(1,2)              &
              .and. j<=indx_couple_grid_center(2,2)              &
              .and. k>=indx_couple_grid_center(1,3)              &
              .and. k<=indx_couple_grid_center(2,3)              &
               ) then
               Vx(i,j,k)= Vx(i,j,k)+ Vx_Far_1D(m)
               Vy(i,j,k)= Vy(i,j,k)+ Vy_Far_1D(m)
               Vz(i,j,k)= Vz(i,j,k)+ Vz_Far_1D(m)
              Txx(i,j,k)=Txx(i,j,k)+Txx_Far_1D(m)
              Tyy(i,j,k)=Tyy(i,j,k)+Tyy_Far_1D(m)
              Tzz(i,j,k)=Tzz(i,j,k)+Tzz_Far_1D(m)
              Tyz(i,j,k)=Tyz(i,j,k)+Tyz_Far_1D(m)
              Txz(i,j,k)=Txz(i,j,k)+Txz_Far_1D(m)
              Txy(i,j,k)=Txy(i,j,k)+Txy_Far_1D(m)
          end if
       end do
       end do
       end do
    end do
!-------------------------------------------------------------------------------
end subroutine add_far2main
!-------------------------------------------------------------------------------

!===============================================================================
subroutine free_charac  !-- characteristic variable free surface boundary
!===============================================================================
    integer i,j,k
    real(SP) v1,v2,v3,t11,t22,t33,t12,t13,t23
    real(SP) lam,miu,lam2mu,rrho,f1,f2,fct
    real(SP) :: Tx,Ty,Tz
    
    if (.not. freenode) return
    
       k=nk2
    do j=nj1,nj2
    do i=ni1,ni2
    
       lam=C13(i,j,k);miu=C66(i,j,k);lam2mu=lam+2.0*miu;
       rrho=rho(i,j,k)
       f1=sqrt(rrho*lam2mu); f2=sqrt(rrho*miu); fct=lam/lam2mu
       v1 = Vx (i,j,k)
       v2 = Vy (i,j,k)
       v3 = Vz (i,j,k)
       t11= Txx(i,j,k)
       t22= Tyy(i,j,k)
       t33= Tzz(i,j,k)
       t12= Txy(i,j,k)
       t13= Txz(i,j,k)
       t23= Tyz(i,j,k)
    
       Tx=-TxSrc(i,j); Ty=-TySrc(i,j); Tz= TzSrc(i,j)
       Vx (i,j,k)=v1-(t13+Tx)/f2
       Vy (i,j,k)=v2-(t23+Ty)/f2
       Vz (i,j,k)=v3-(t33-Tz)/f1
       Txx(i,j,k)=t11-fct*(t33-Tz)
       Tyy(i,j,k)=t22-fct*(t33-Tz)
       Tzz(i,j,k)=Tz
       Txz(i,j,k)=Tx
       Tyz(i,j,k)=Ty
    end do
    end do
!===============================================================================
end subroutine free_charac
!===============================================================================

!===============================================================================
subroutine free_extrap  !--- extrapolate wavefield on grids above surface
!===============================================================================
    integer i,j,k
    
    if (freenode) then
    do k=nk2+1,nz2
    do j=nj1,nj2
    do i=ni1,ni2
       Vx (i,j,k)=4.0*Vx (i,j,k-1)-6.0*Vx (i,j,k-2)+4.0*Vx (i,j,k-3)-Vx (i,j,k-4)
       Vy (i,j,k)=4.0*Vy (i,j,k-1)-6.0*Vy (i,j,k-2)+4.0*Vy (i,j,k-3)-Vy (i,j,k-4)
       Vz (i,j,k)=4.0*Vz (i,j,k-1)-6.0*Vz (i,j,k-2)+4.0*Vz (i,j,k-3)-Vz (i,j,k-4)
       Txx(i,j,k)=4.0*Txx(i,j,k-1)-6.0*Txx(i,j,k-2)+4.0*Txx(i,j,k-3)-Txx(i,j,k-4)
       Tzz(i,j,k)=4.0*Tzz(i,j,k-1)-6.0*Tzz(i,j,k-2)+4.0*Tzz(i,j,k-3)-Tzz(i,j,k-4)
       Tyy(i,j,k)=4.0*Tyy(i,j,k-1)-6.0*Tyy(i,j,k-2)+4.0*Tyy(i,j,k-3)-Tyy(i,j,k-4)
       Txy(i,j,k)=4.0*Txy(i,j,k-1)-6.0*Txy(i,j,k-2)+4.0*Txy(i,j,k-3)-Txy(i,j,k-4)
       Txz(i,j,k)=4.0*Txz(i,j,k-1)-6.0*Txz(i,j,k-2)+4.0*Txz(i,j,k-3)-Txz(i,j,k-4)
       Tyz(i,j,k)=4.0*Tyz(i,j,k-1)-6.0*Tyz(i,j,k-2)+4.0*Tyz(i,j,k-3)-Tyz(i,j,k-4)
    end do
    end do
    end do
    end if
!===============================================================================
end subroutine free_extrap
!===============================================================================

!===============================================================================
subroutine free_vext  !--- only extrapolate velocity above surface
!===============================================================================
    integer i,j,k
    
    if (freenode) then
    do k=nk2+1,nz2
    do j=nj1,nj2
    do i=ni1,ni2
       Vx (i,j,k)=4.0*Vx (i,j,k-1)-6.0*Vx (i,j,k-2)+4.0*Vx (i,j,k-3)-Vx (i,j,k-4)
       Vy (i,j,k)=4.0*Vy (i,j,k-1)-6.0*Vy (i,j,k-2)+4.0*Vy (i,j,k-3)-Vy (i,j,k-4)
       Vz (i,j,k)=4.0*Vz (i,j,k-1)-6.0*Vz (i,j,k-2)+4.0*Vz (i,j,k-3)-Vz (i,j,k-4)
    end do
    end do
    end do
    end if
!===============================================================================
end subroutine free_vext
!===============================================================================

!===============================================================================
subroutine free_timg  !--- stress imaging of stress
!===============================================================================
    integer i,j,k
    
    if (freenode) then
    do k=1,SEIS_GEO
    do j=nj1,nj2
    do i=ni1,ni2
       Tzz(i,j,nk2+k)=2.0_SP*TzSrc(i,j)-Tzz(i,j,nk2-k)
       Txz(i,j,nk2+k)=2.0_SP*TxSrc(i,j)-Txz(i,j,nk2-k)
       Tyz(i,j,nk2+k)=2.0_SP*TySrc(i,j)-Tyz(i,j,nk2-k)
    end do
    end do
    end do
    
    do j=nj1,nj2
    do i=ni1,ni2
       Tzz(i,j,nk2)=TzSrc(i,j)
       Txz(i,j,nk2)=TxSrc(i,j)
       Tyz(i,j,nk2)=TySrc(i,j)
    end do
    end do
    end if
!===============================================================================
end subroutine free_timg
!===============================================================================

!===============================================================================
subroutine macdrp_mesg_init  !--- initial the continuous communication
!===============================================================================
    call mesg_init_LxF
    call mesg_init_LxB
    call mesg_init_LyF
    call mesg_init_LyB
    call mesg_init_LzF
    call mesg_init_LzB
!===============================================================================
end subroutine macdrp_mesg_init
!===============================================================================

!===============================================================================
!--- message for each direction and forward or backward operator
!===============================================================================

!===============================================================================
subroutine mesg_init_LxF
!===============================================================================
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufX1,NBufXL,SEISMPI_DATATYPE,neigid(1,1),1131,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevX2,NBufXL,SEISMPI_DATATYPE,neigid(1,2),1131,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufX2,NBufXS,SEISMPI_DATATYPE,neigid(1,2),1211,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevX1,NBufXS,SEISMPI_DATATYPE,neigid(1,1),1211,SWMPI_COMM,r2,ierr)
reqXF=(/ s1,r1,s2,r2 /)
#else
! to X1
call MPI_SEND_INIT(Txx(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1131,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1132,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1135,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni1,ny1,nz1),1,DTypeXL,neigid(1,1),1139,SWMPI_COMM,s9,ierr)
! from X2
call MPI_RECV_INIT(Txx(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1131,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1132,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1135,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni2+1,ny1,nz1),1,DTypeXL,neigid(1,2),1139,SWMPI_COMM,r9,ierr)
! put into array
!reqXF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXF(1:12)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)

!to X2
call MPI_SEND_INIT(Txx(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1211,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1212,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1215,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-LenFDS+1,ny1,nz1),1,DTypeXS,neigid(1,2),1219,SWMPI_COMM,s9,ierr)
!from X1
call MPI_RECV_INIT(Txx(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1211,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1212,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1215,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-LenFDS,ny1,nz1),1,DTypeXS,neigid(1,1),1219,SWMPI_COMM,r9,ierr)
! put into array
!reqXF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXF(13:24)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)
#endif
!===============================================================================
end subroutine mesg_init_LxF
!===============================================================================

!===============================================================================
subroutine mesg_init_LxB
!===============================================================================
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufX1,NBufXS,SEISMPI_DATATYPE,neigid(1,1),1111,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevX2,NBufXS,SEISMPI_DATATYPE,neigid(1,2),1111,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufX2,NBufXL,SEISMPI_DATATYPE,neigid(1,2),1231,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevX1,NBufXL,SEISMPI_DATATYPE,neigid(1,1),1231,SWMPI_COMM,r2,ierr)
reqXB=(/ s1,r1,s2,r2 /)
#else
! to X1
call MPI_SEND_INIT(Txx(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1111,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1112,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1115,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni1,ny1,nz1),1,DTypeXS,neigid(1,1),1119,SWMPI_COMM,s9,ierr)
! from X2
call MPI_RECV_INIT(Txx(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1111,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1112,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1115,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni2+1,ny1,nz1),1,DTypeXS,neigid(1,2),1119,SWMPI_COMM,r9,ierr)
! put into array
!reqXB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXB(1:12)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)

! to X2
call MPI_SEND_INIT(Txx(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1231,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1232,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1235,SWMPI_COMM,s5,ierr)
!call MPI_SEND_INIT(Tyz(ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (ni2-LenFDL+1,ny1,nz1),1,DTypeXL,neigid(1,2),1239,SWMPI_COMM,s9,ierr)
! from  X1
call MPI_RECV_INIT(Txx(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1231,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1232,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1235,SWMPI_COMM,r5,ierr)
!call MPI_RECV_INIT(Tyz(ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (ni1-LenFDL,ny1,nz1),1,DTypeXL,neigid(1,1),1239,SWMPI_COMM,r9,ierr)
! put into array
!reqXB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqXB(13:24)=(/ s1,r1,s4,r4,s5,r5,s7,r7,s8,r8,s9,r9 /)
#endif
!===============================================================================
end subroutine mesg_init_LxB
!===============================================================================

!===============================================================================
subroutine mesg_init_LyF
!===============================================================================
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufY1,NBufYL,SEISMPI_DATATYPE,neigid(2,1),2131,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevY2,NBufYL,SEISMPI_DATATYPE,neigid(2,2),2131,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufY2,NBufYS,SEISMPI_DATATYPE,neigid(2,2),2211,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevY1,NBufYS,SEISMPI_DATATYPE,neigid(2,1),2211,SWMPI_COMM,r2,ierr)
reqYF=(/ s1,r1,s2,r2 /)
#else
! to Y1
!call MPI_SEND_INIT(Txx(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2131,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2132,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2133,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2134,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj1,nz1),1,DTypeYL,neigid(2,1),2139,SWMPI_COMM,s9,ierr)
! from Y2
!call MPI_RECV_INIT(Txx(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2131,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2132,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2133,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2134,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj2+1,nz1),1,DTypeYL,neigid(2,2),2139,SWMPI_COMM,r9,ierr)
! put into array
!reqYF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYF(1:12)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
!-----

! to Y2
!call MPI_SEND_INIT(Txx(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2211,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2212,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2213,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2214,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-LenFDS+1,nz1),1,DTypeYS,neigid(2,2),2219,SWMPI_COMM,s9,ierr)
! from Y1
!call MPI_RECV_INIT(Txx(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2211,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2212,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2213,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2214,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-LenFDS,nz1),1,DTypeYS,neigid(2,1),2219,SWMPI_COMM,r9,ierr)
! put into array
!reqYF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYF(13:24)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
#endif
!===============================================================================
end subroutine mesg_init_LyF
!===============================================================================

!===============================================================================
subroutine mesg_init_LyB
!===============================================================================
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufY1,NBufYS,SEISMPI_DATATYPE,neigid(2,1),2111,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevY2,NBufYS,SEISMPI_DATATYPE,neigid(2,2),2111,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufY2,NBufYL,SEISMPI_DATATYPE,neigid(2,2),2231,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevY1,NBufYL,SEISMPI_DATATYPE,neigid(2,1),2231,SWMPI_COMM,r2,ierr)
reqYB=(/ s1,r1,s2,r2 /)
#else
! to Y1
!call MPI_SEND_INIT(Txx(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2111,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2112,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2113,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2114,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj1,nz1),1,DTypeYS,neigid(2,1),2119,SWMPI_COMM,s9,ierr)
! from Y2
!call MPI_RECV_INIT(Txx(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2111,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2112,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2113,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2114,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj2+1,nz1),1,DTypeYS,neigid(2,2),2119,SWMPI_COMM,r9,ierr)
! put into array
!reqYB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYB(1:12)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)

! to Y2
!call MPI_SEND_INIT(Txx(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2231,SWMPI_COMM,s1,ierr)
call MPI_SEND_INIT(Tyy(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2232,SWMPI_COMM,s2,ierr)
!call MPI_SEND_INIT(Tzz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2233,SWMPI_COMM,s3,ierr)
call MPI_SEND_INIT(Txy(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2234,SWMPI_COMM,s4,ierr)
!call MPI_SEND_INIT(Txz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,nj2-LenFDL+1,nz1),1,DTypeYL,neigid(2,2),2239,SWMPI_COMM,s9,ierr)
! from Y1
!call MPI_RECV_INIT(Txx(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2231,SWMPI_COMM,r1,ierr)
call MPI_RECV_INIT(Tyy(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2232,SWMPI_COMM,r2,ierr)
!call MPI_RECV_INIT(Tzz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2233,SWMPI_COMM,r3,ierr)
call MPI_RECV_INIT(Txy(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2234,SWMPI_COMM,r4,ierr)
!call MPI_RECV_INIT(Txz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,nj1-LenFDL,nz1),1,DTypeYL,neigid(2,1),2239,SWMPI_COMM,r9,ierr)
! put into array
!reqYB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqYB(13:24)=(/ s2,r2,s4,r4,s6,r6,s7,r7,s8,r8,s9,r9 /)
#endif
!===============================================================================
end subroutine mesg_init_LyB
!===============================================================================

!===============================================================================
subroutine mesg_init_LzF
!===============================================================================
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufZ1,NBufZL,SEISMPI_DATATYPE,neigid(3,1),3131,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevZ2,NBufZL,SEISMPI_DATATYPE,neigid(3,2),3131,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufZ2,NBufZS,SEISMPI_DATATYPE,neigid(3,2),3211,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevZ1,NBufZS,SEISMPI_DATATYPE,neigid(3,1),3211,SWMPI_COMM,r2,ierr)
reqZF=(/ s1,r1,s2,r2 /)
#else
!to Z1
!call MPI_SEND_INIT(Txx(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3131,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3132,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3133,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3134,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3135,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3136,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3137,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3138,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk1),1,DTypeZL,neigid(3,1),3139,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3131,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3132,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3133,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3134,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3135,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3136,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3137,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3138,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk2+1),1,DTypeZL,neigid(3,2),3139,SWMPI_COMM,r9,ierr)
! put into array
!reqZF(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZF(1:12)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to Z2
!call MPI_SEND_INIT(Txx(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3211,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3212,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3213,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3214,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3215,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3216,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3217,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3218,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-LenFDS+1),1,DTypeZS,neigid(3,2),3219,SWMPI_COMM,s9,ierr)
!from Z1
!call MPI_RECV_INIT(Txx(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3211,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3212,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3213,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3214,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3215,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3216,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3217,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3218,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-LenFDS),1,DTypeZS,neigid(3,1),3219,SWMPI_COMM,r9,ierr)
! put into array
!reqZF(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZF(13:24)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
#endif
!===============================================================================
end subroutine mesg_init_LzF
!===============================================================================

!===============================================================================
subroutine mesg_init_LzB
!===============================================================================
integer s1,s2,s3,s4,s5,s6,s7,s8,s9
integer r1,r2,r3,r4,r5,r6,r7,r8,r9
#ifdef MPIBuffered
call MPI_SEND_INIT(BufZ1,NBufZS,SEISMPI_DATATYPE,neigid(3,1),3111,SWMPI_COMM,s1,ierr)
call MPI_RECV_INIT(RevZ2,NBufZS,SEISMPI_DATATYPE,neigid(3,2),3111,SWMPI_COMM,r1,ierr)
call MPI_SEND_INIT(BufZ2,NBufZL,SEISMPI_DATATYPE,neigid(3,2),3231,SWMPI_COMM,s2,ierr)
call MPI_RECV_INIT(RevZ1,NBufZL,SEISMPI_DATATYPE,neigid(3,1),3231,SWMPI_COMM,r2,ierr)
reqZB=(/ s1,r1,s2,r2 /)
#else
!to Z1
!call MPI_SEND_INIT(Txx(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3111,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3112,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3113,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3114,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3115,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3116,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3117,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3118,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk1),1,DTypeZS,neigid(3,1),3119,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3111,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3112,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3113,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3114,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3115,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3116,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3117,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3118,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk2+1),1,DTypeZS,neigid(3,2),3119,SWMPI_COMM,r9,ierr)
! put into array
!reqZB(1:18)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZB(1:12)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)

!to Z2
!call MPI_SEND_INIT(Txx(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3231,SWMPI_COMM,s1,ierr)
!call MPI_SEND_INIT(Tyy(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3232,SWMPI_COMM,s2,ierr)
call MPI_SEND_INIT(Tzz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3233,SWMPI_COMM,s3,ierr)
!call MPI_SEND_INIT(Txy(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3234,SWMPI_COMM,s4,ierr)
call MPI_SEND_INIT(Txz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3235,SWMPI_COMM,s5,ierr)
call MPI_SEND_INIT(Tyz(nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3236,SWMPI_COMM,s6,ierr)
call MPI_SEND_INIT(Vx (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3237,SWMPI_COMM,s7,ierr)
call MPI_SEND_INIT(Vy (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3238,SWMPI_COMM,s8,ierr)
call MPI_SEND_INIT(Vz (nx1,ny1,nk2-LenFDL+1),1,DTypeZL,neigid(3,2),3239,SWMPI_COMM,s9,ierr)
!from Z2
!call MPI_RECV_INIT(Txx(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3231,SWMPI_COMM,r1,ierr)
!call MPI_RECV_INIT(Tyy(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3232,SWMPI_COMM,r2,ierr)
call MPI_RECV_INIT(Tzz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3233,SWMPI_COMM,r3,ierr)
!call MPI_RECV_INIT(Txy(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3234,SWMPI_COMM,r4,ierr)
call MPI_RECV_INIT(Txz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3235,SWMPI_COMM,r5,ierr)
call MPI_RECV_INIT(Tyz(nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3236,SWMPI_COMM,r6,ierr)
call MPI_RECV_INIT(Vx (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3237,SWMPI_COMM,r7,ierr)
call MPI_RECV_INIT(Vy (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3238,SWMPI_COMM,r8,ierr)
call MPI_RECV_INIT(Vz (nx1,ny1,nk1-LenFDL),1,DTypeZL,neigid(3,1),3239,SWMPI_COMM,r9,ierr)
! put into array
!reqZB(19:36)=(/ s1,r1,s2,r2,s3,r3,s4,r4,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
reqZB(13:24)=(/ s3,r3,s5,r5,s6,r6,s7,r7,s8,r8,s9,r9 /)
#endif
!===============================================================================
end subroutine mesg_init_LzB
!===============================================================================

#ifdef MPIBuffered
!-------------------------------------------------------------------------------
!--- buffered communication
!-------------------------------------------------------------------------------
subroutine fill_buff_LxF
integer :: i,j,k,n
! to X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX1(n)=Tyy(i,j,k)
   !n=n+1; BufX1(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDL-1
   n=n+1; BufX1(n)=Txz(i,j,k)
   !n=n+1; BufX1(n)=Tyz(i,j,k)
end do
end do
end do
! to X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX2(n)=Tyy(i,j,k)
   !n=n+1; BufX2(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDS+1,ni2
   n=n+1; BufX2(n)=Txz(i,j,k)
   !n=n+1; BufX2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LxF
!-------------------------------------------------------------------------------
subroutine recv_buff_LxF
integer :: i,j,k,n
! from X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1;  Vx(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1;  Vy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1;  Vz(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1; Txx(i,j,k)=RevX1(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX1(n)
   !n=n+1; Tzz(i,j,k)=RevX1(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1; Txy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDS,ni1-1
   n=n+1; Txz(i,j,k)=RevX1(n)
   !n=n+1; Tyz(i,j,k)=RevX1(n)
end do
end do
end do
! from X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1;  Vx(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1;  Vy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1;  Vz(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1; Txx(i,j,k)=RevX2(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX2(n)
   !n=n+1; Tzz(i,j,k)=RevX2(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1; Txy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDL
   n=n+1; Txz(i,j,k)=RevX2(n)
   !n=n+1; Tyz(i,j,k)=RevX2(n)
end do
end do
end do
end subroutine recv_buff_LxF
!-------------------------------------------------------------------------------
subroutine fill_buff_LxB
integer :: i,j,k,n
! to X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX1(n)=Tyy(i,j,k)
   !n=n+1; BufX1(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1,ni1+LenFDS-1
   n=n+1; BufX1(n)=Txz(i,j,k)
   !n=n+1; BufX1(n)=Tyz(i,j,k)
end do
end do
end do
! to X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)=Txx(i,j,k)
end do
end do
end do
   !n=n+1; BufX2(n)=Tyy(i,j,k)
   !n=n+1; BufX2(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2-LenFDL+1,ni2
   n=n+1; BufX2(n)=Txz(i,j,k)
   !n=n+1; BufX2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LxB
!-------------------------------------------------------------------------------
subroutine recv_buff_LxB
integer :: i,j,k,n
! from X1
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1;  Vx(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1;  Vy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1;  Vz(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1; Txx(i,j,k)=RevX1(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX1(n)
   !n=n+1; Tzz(i,j,k)=RevX1(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1; Txy(i,j,k)=RevX1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni1-LenFDL,ni1-1
   n=n+1; Txz(i,j,k)=RevX1(n)
   !n=n+1; Tyz(i,j,k)=RevX1(n)
end do
end do
end do
! from X2
n=0
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1;  Vx(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1;  Vy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1;  Vz(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1; Txx(i,j,k)=RevX2(n)
end do
end do
end do
   !n=n+1; Tyy(i,j,k)=RevX2(n)
   !n=n+1; Tzz(i,j,k)=RevX2(n)
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1; Txy(i,j,k)=RevX2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj2
do i=ni2+1,ni2+LenFDS
   n=n+1; Txz(i,j,k)=RevX2(n)
   !n=n+1; Tyz(i,j,k)=RevX2(n)
end do
end do
end do
end subroutine recv_buff_LxB
!-------------------------------------------------------------------------------
subroutine fill_buff_LyF
integer :: i,j,k,n
! to Y1
n=0
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vz(i,j,k)
end do
end do
end do
   !n=n+1; BufY1(n)=Txx(i,j,k)
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)=Tyy(i,j,k)
end do
end do
end do
   !n=n+1; BufY1(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)=Txy(i,j,k)
end do
end do
end do
   !n=n+1; BufY1(n)=Txz(i,j,k)
do k=nk1,nk2
do j=nj1,nj1+LenFDL-1
do i=ni1,ni2
   n=n+1; BufY1(n)=Tyz(i,j,k)
end do
end do
end do
! to Y2
n=0
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vz(i,j,k)
end do
end do
end do
   !n=n+1; BufY2(n)=Txx(i,j,k)
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)=Tyy(i,j,k)
end do
end do
end do
   !n=n+1; BufY2(n)=Tzz(i,j,k)
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)=Txy(i,j,k)
end do
end do
end do
   !n=n+1; BufY2(n)=Txz(i,j,k)
do k=nk1,nk2
do j=nj2-LenFDS+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LyF
!-------------------------------------------------------------------------------
subroutine recv_buff_LyF
integer :: i,j,k,n
! from Y1
n=0
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY1(n)
   n=n+1; Tyy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY1(n)
   n=n+1; Txy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDS,nj1-1
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY1(n)
   n=n+1; Tyz(i,j,k)=RevY1(n)
end do
end do
end do
! from Y2
n=0
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY2(n)
   n=n+1; Tyy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY2(n)
   n=n+1; Txy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDL
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY2(n)
   n=n+1; Tyz(i,j,k)=RevY2(n)
end do
end do
end do
end subroutine recv_buff_LyF
!-------------------------------------------------------------------------------
subroutine fill_buff_LyB
integer :: i,j,k,n
! to Y1
n=0
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   n=n+1; BufY1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   !n=n+1; BufY1(n)=Txx(i,j,k)
   n=n+1; BufY1(n)=Tyy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   !n=n+1; BufY1(n)=Tzz(i,j,k)
   n=n+1; BufY1(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj1,nj1+LenFDS-1
do i=ni1,ni2
   !n=n+1; BufY1(n)=Txz(i,j,k)
   n=n+1; BufY1(n)=Tyz(i,j,k)
end do
end do
end do
! to Y2
n=0
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   n=n+1; BufY2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   !n=n+1; BufY2(n)=Txx(i,j,k)
   n=n+1; BufY2(n)=Tyy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   !n=n+1; BufY2(n)=Tzz(i,j,k)
   n=n+1; BufY2(n)=Txy(i,j,k)
end do
end do
end do
do k=nk1,nk2
do j=nj2-LenFDL+1,nj2
do i=ni1,ni2
   !n=n+1; BufY2(n)=Txz(i,j,k)
   n=n+1; BufY2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LyB
!-------------------------------------------------------------------------------
subroutine recv_buff_LyB
integer :: i,j,k,n
! from Y1
n=0
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY1(n)
   n=n+1; Tyy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY1(n)
   n=n+1; Txy(i,j,k)=RevY1(n)
end do
end do
end do
do k=nk1,nk2
do j=nj1-LenFDL,nj1-1
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY1(n)
   n=n+1; Tyz(i,j,k)=RevY1(n)
end do
end do
end do
! from Y2
n=0
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevY2(n)
   n=n+1; Tyy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   !n=n+1; Tzz(i,j,k)=RevY2(n)
   n=n+1; Txy(i,j,k)=RevY2(n)
end do
end do
end do
do k=nk1,nk2
do j=nj2+1,nj2+LenFDS
do i=ni1,ni2
   !n=n+1; Txz(i,j,k)=RevY2(n)
   n=n+1; Tyz(i,j,k)=RevY2(n)
end do
end do
end do
end subroutine recv_buff_LyB
!-------------------------------------------------------------------------------
subroutine fill_buff_LzF
integer :: i,j,k,n
! to Z1
n=0
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txx(i,j,k)
   !n=n+1; BufZ1(n)=Tyy(i,j,k)
   n=n+1; BufZ1(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txy(i,j,k)
   n=n+1; BufZ1(n)=Txz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDL-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)=Tyz(i,j,k)
end do
end do
end do
! to Z2
n=0
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txx(i,j,k)
   !n=n+1; BufZ2(n)=Tyy(i,j,k)
   n=n+1; BufZ2(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txy(i,j,k)
   n=n+1; BufZ2(n)=Txz(i,j,k)
end do
end do
end do
do k=nk2-LenFDS+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LzF
!-------------------------------------------------------------------------------
subroutine recv_buff_LzF
integer :: i,j,k,n
! from Z1
n=0
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ1(n)
   !n=n+1; Tyy(i,j,k)=RevZ1(n)
   n=n+1; Tzz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ1(n)
   n=n+1; Txz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDS,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ1(n)
end do
end do
end do
! from Z2
n=0
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ2(n)
   !n=n+1; Tyy(i,j,k)=RevZ2(n)
   n=n+1; Tzz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ2(n)
   n=n+1; Txz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDL
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ2(n)
end do
end do
end do
end subroutine recv_buff_LzF
!-------------------------------------------------------------------------------
subroutine fill_buff_LzB
integer :: i,j,k,n
! to Z1
n=0
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vx(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vy(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)= Vz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txx(i,j,k)
   !n=n+1; BufZ1(n)=Tyy(i,j,k)
   n=n+1; BufZ1(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ1(n)=Txy(i,j,k)
   n=n+1; BufZ1(n)=Txz(i,j,k)
end do
end do
end do
do k=nk1,nk1+LenFDS-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ1(n)=Tyz(i,j,k)
end do
end do
end do
! to Z2
n=0
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vx(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vy(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)= Vz(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txx(i,j,k)
   !n=n+1; BufZ2(n)=Tyy(i,j,k)
   n=n+1; BufZ2(n)=Tzz(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; BufZ2(n)=Txy(i,j,k)
   n=n+1; BufZ2(n)=Txz(i,j,k)
end do
end do
end do
do k=nk2-LenFDL+1,nk2
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; BufZ2(n)=Tyz(i,j,k)
end do
end do
end do
end subroutine fill_buff_LzB
!-------------------------------------------------------------------------------
subroutine recv_buff_LzB
integer :: i,j,k,n
! from Z1
n=0
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ1(n)
   !n=n+1; Tyy(i,j,k)=RevZ1(n)
   n=n+1; Tzz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ1(n)
   n=n+1; Txz(i,j,k)=RevZ1(n)
end do
end do
end do
do k=nk1-LenFDL,nk1-1
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ1(n)
end do
end do
end do
! from Z2
n=0
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vx(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vy(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1;  Vz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txx(i,j,k)=RevZ2(n)
   !n=n+1; Tyy(i,j,k)=RevZ2(n)
   n=n+1; Tzz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   !n=n+1; Txy(i,j,k)=RevZ2(n)
   n=n+1; Txz(i,j,k)=RevZ2(n)
end do
end do
end do
do k=nk2+1,nk2+LenFDS
do j=nj1,nj2
do i=ni1,ni2
   n=n+1; Tyz(i,j,k)=RevZ2(n)
end do
end do
end do
end subroutine recv_buff_LzB
!-------------------------------------------------------------------------------
#endif

!-------------------------------------------------------------------------------
function macdrp_get_fnm_coupling(pnm,n_i,n_j,n_k) result(filenm)
!-------------------------------------------------------------------------------
    character (len=*) :: pnm
    integer,intent(in) :: n_i,n_j,n_k
    character (len=SEIS_STRLEN) :: filenm
    filenm=trim(pnm)//'/'//'coupling'  &
           //'_'//trim(set_mpi_subfix(n_i,n_j,n_k))//'.bin'
!-------------------------------------------------------------------------------
end function macdrp_get_fnm_coupling
!-------------------------------------------------------------------------------

!===============================================================================
subroutine macdrp_print_info(fid) 
!===============================================================================

!-------------------------------------------------------------------------------
!-- Description:
!-------------------------------------------------------------------------------
!--   print variables and compiler macro

!-------------------------------------------------------------------------------
!-- Input arguments:
!-------------------------------------------------------------------------------
!--
    integer,intent(in) :: fid        !-- id of opened info file

!-------------------------------------------------------------------------------
!-- Local variables:
!-------------------------------------------------------------------------------

    character (len=300) :: FFLAG     !-- compiler flag

!-------------------------------------------------------------------------------
!-- Entry point: 
!-------------------------------------------------------------------------------

    write(fid,*)
    write(fid,"(a)") "#========================================================"
    write(fid,"(a)") "Compiler MACRO used for macdrp_mod:"
    write(fid,"(a)") "#========================================================"
    FFLAG = ''

!---- MPI or feather related
#ifdef MPIBuffered
    FFLAG = trim(FFLAG)//' -DMPIBuffered'
#endif
#ifdef MPIBARRIER
    FFLAG = trim(FFLAG)//' -DMPIBARRIER'
#endif
#ifdef VERBOSE
    FFLAG = trim(FFLAG)//' -DVERBOSE'
#endif

!---- medium related
#ifdef AnisGene
    FFLAG = trim(FFLAG)//' -DAnisGene'
#endif
#ifdef AnisVTI
    FFLAG = trim(FFLAG)//' -DAnisVTI'
#endif
#ifdef WATER
    FFLAG = trim(FFLAG)//' -DWATER'
#endif

!---- free surface boundary related
#ifdef CondFreeCharac
    FFLAG = trim(FFLAG)//' -DCondFreeCharac'
#endif
#ifdef CondFreeTIMG
    FFLAG = trim(FFLAG)//' -DCondFreeTIMG'
#endif
#ifdef CondFreeVEXT
    FFLAG = trim(FFLAG)//' -DCondFreeVEXT'
#endif
#ifdef CondFreeVHOC
    FFLAG = trim(FFLAG)//' -DCondFreeVHOC'
#endif
#ifdef CondFreeVLOW
    FFLAG = trim(FFLAG)//' -DCondFreeVLOW'
#endif

    !---- write FFLAG to the output file -----
    write(fid,"(a)") "FFLAG = "//trim(FFLAG)
    write(fid,*)

    write(fid,"(a)") "#========================================================"
    write(fid,"(a)") "Input parameter values in macdrp_mod:"
    write(fid,"(a)") "#========================================================"
    write(fid,"(a)") "  no input parameter"
    write(fid,*)

!===============================================================================
end subroutine macdrp_print_info
!===============================================================================

end module macdrp_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
