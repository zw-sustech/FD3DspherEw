divert(-1)dnl
dnl macro to generate mod_abs_npml.F90
changequote([,])dnl
divert(0)dnl

program seis3d_wave

!-------------------------------------------------------------------------------
! Description:
!-------------------------------------------------------------------------------
!
! This is the main program to simulate seismic wave propagation
!
! Author: Wei ZHANG     Email: zhangwei.zw@gmail.com
! Copyright (C) 2006 Wei ZHANG
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
! Dependent modules:
!-------------------------------------------------------------------------------
    use mpi
    use para_mod
    use io_mod
    use mpi_mod
    use grid_mod
    use media_mod
    use src_mod
    use abs_mod
    use macdrp_mod
#ifdef WithOMP
    use omp_lib
#endif

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
    implicit none
    integer :: ntime,ierr
    integer :: fid_info  !-- fid for info file

!-------------------------------------------------------------------------------
! Program Entries:
!-------------------------------------------------------------------------------

    call MPI_INIT(ierr)
    
    call get_conf_name(fnm_conf)
    
    call swmpi_init(fnm_conf)
    call swmpi_cart_creat
    
    call para_init(fnm_conf)
    call swmpi_reinit_para
    
    call grid_fnm_init(fnm_conf)
    call grid_alloc
    call grid_coord_import(thisid(1),thisid(2),thisid(3))
    call grid_metric_import(thisid(1),thisid(2),thisid(3))
    
    call media_fnm_init(fnm_conf)
    call media_alloc
    call media_import(thisid(1),thisid(2),thisid(3))
    
    call src_fnm_init(fnm_conf)
    call src_import(thisid(1),thisid(2),thisid(3))
    !call src_choose
    
    call io_init(fnm_conf)
    call io_snap_read(fnm_conf)
    call io_snap_locate(thisid(1),thisid(2),thisid(3))
    call io_pt_import
    call io_seismo_init
    
    !call grid_dealloc(iscoord=.true.)
    
    call macdrp_init(fnm_conf)
    
    call abs_init(fnm_conf)
    
    call swmpi_datatype
    call macdrp_mesg_init
    
    ntime=0
    
    call io_rest_import(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
    
    call swmpi_time_init(filename_log,ntime)

#ifdef WithOMP
    call OMP_SET_NUM_THREADS(2)
#endif

!-------------------------------------------------------------------------------
! print info in each module
!-------------------------------------------------------------------------------
    if ( masternode ) then
    
        fid_info = 9101
    
        !-- open info file --
        call io_info_file_open(myid,fid_info)
    
        !-- write info --
        call io_print_info    (fid_info)
        call macdrp_print_info(fid_info)
    
        !-- close info file
        call io_info_file_close(fid_info)
    end if

!-------------------------------------------------------------------------------
! seis3d_wave_8op_p3
!-------------------------------------------------------------------------------
! 5-1B BBB
! 8-4B FFB
! 1-1A FFF
! 4-4A BBF

! 7-3B BFB
! 6-2B FBB
! 3-3A FBF
! 2-2A BFF

!===============================================================================
loop_time: do
!===============================================================================

    if ( ntime>nt ) exit

divert([-1])
define([_UPDATEONESTEP],
[dnl
!-------------------------------------------------------------------------------
!-- [$3]: [$1] - [$2]
!-------------------------------------------------------------------------------
!-- prepare
    call swmpi_time_write(ntime,filename_log)
    call macdrp_syn
    call abs_syn
!-- the 1th stage
    call set_cur_time(ntime,0.0_SP)
#ifdef SrcSurface
    call src_surface(ntime,0.0_SP,stept)
#endif
    call macdrp_[$1]
    call abs_[$1]
      call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,0.0_SP,stept)
      call src_force(hVx,hVy,hVz,ntime,0.0_SP,stept)
    call macdrp_RK_beg(firRKa(1),firRKb(1))
    call abs_RK_beg(firRKa(1),firRKb(1))
!-- the 2th stage
    call set_cur_time(ntime,firRKa(1))
#ifdef SrcSurface
    call src_surface(ntime,firRKa(1),stept)
#endif
    call macdrp_[$2]
    call abs_[$2]
      call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(1),stept)
      call src_force(hVx,hVy,hVz,ntime,firRKa(1),stept)
    call macdrp_RK_inn(firRKa(2),firRKb(2))
    call abs_RK_inn(firRKa(2),firRKb(2))
!-- the 3th stage
    call set_cur_time(ntime,firRKa(2))
#ifdef SrcSurface
    call src_surface(ntime,firRKa(2),stept)
#endif
    call macdrp_[$1]
    call abs_[$1]
      call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(2),stept)
      call src_force(hVx,hVy,hVz,ntime,firRKa(2),stept)
    call macdrp_RK_inn(firRKa(3),firRKb(3))
    call abs_RK_inn(firRKa(3),firRKb(3))
!-- the 4th stage
    call set_cur_time(ntime,firRKa(3))
#ifdef SrcSurface
    call src_surface(ntime,firRKa(3),stept)
#endif
    call macdrp_[$2]
    call abs_[$2]
      call src_stress(hTxx,hTyy,hTzz,hTxy,hTxz,hTyz,ntime,firRKa(3),stept)
      call src_force(hVx,hVy,hVz,ntime,firRKa(3),stept)
    call macdrp_RK_fin(firRKb(4))
    call abs_RK_fin(firRKb(4))
#ifdef WITHQS
    call atten_graves
#endif

!-- save result
    ntime=ntime+1
    call macdrp_check(ntime)
    call io_seismo_put(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime)
    call io_wave_export(Vx,Vy,Vz,Txx,Tyy,Tzz,Txy,Txz,Tyz,ntime,stept)
    call io_rest_export(Txx,Tyy,Tzz,Txy,Txz,Tyz,Vx,Vy,Vz,ntime)
])dnl
divert([0])dnl

_UPDATEONESTEP(LxB_LyB_LzB,LxF_LyF_LzF,[1])

_UPDATEONESTEP(LxF_LyF_LzB,LxB_LyB_LzF,[2])

_UPDATEONESTEP(LxF_LyF_LzF,LxB_LyB_LzB,[3])

_UPDATEONESTEP(LxB_LyB_LzF,LxF_LyF_LzB,[4])

_UPDATEONESTEP(LxB_LyF_LzB,LxF_LyB_LzF,[5])

_UPDATEONESTEP(LxF_LyB_LzB,LxB_LyF_LzF,[6])

_UPDATEONESTEP(LxF_LyB_LzF,LxB_LyF_LzB,[7])

_UPDATEONESTEP(LxB_LyF_LzF,LxF_LyB_LzB,[8])

!===============================================================================
end do loop_time
!===============================================================================

    call io_seismo_close
    call io_wave_close
    call swmpi_time_end(filename_log)
    
    call macdrp_destroy
    call grid_dealloc
    call media_destroy
    !call abs_destroy
    call src_destroy
    
    call MPI_FINALIZE(ierr)

!-------------------------------------------------------------------------------
!contains
!-------------------------------------------------------------------------------

end program seis3d_wave

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
! vim:ft=m4:ts=4:sw=4:nu:et:ai:
