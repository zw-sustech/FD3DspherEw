program tool_coupling

!-------------------------------------------------------------------------------
! Description:
!-------------------------------------------------------------------------------
!
! This program inits and distributes the seismic source
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
!#include "mod_macdrp.h"

!-------------------------------------------------------------------------------
! Dependent modules:
!-------------------------------------------------------------------------------
    use constants_mod
    use string_mod
    use nfseis_mod
    use para_mod
    use mpi_mod
    use macdrp_mod, only : macdrp_get_fnm_coupling

!-------------------------------------------------------------------------------
! local variables:
!-------------------------------------------------------------------------------
    implicit none

    logical                          :: &
        gflag_x1,gflag_x2,gflag_y1,gflag_y2,gflag_z1,gflag_z2 !- if coupling at the face
    integer,dimension(2,SEIS_GEO)    :: &
        gindx_x1,gindx_x2,gindx_y1,gindx_y2,gindx_z1,gindx_z2 !- index of coupling grid
    character (len=SEIS_STRLEN)      :: &
        gfnm_x1,gfnm_x2,gfnm_y1,gfnm_y2,gfnm_z1,gfnm_z2       !- far field file
    integer ::  &
        ntwin_Far_field,             &  !- total steps of far field simulation (double local steps)
        max_thread_couple, &  !- maximum coupling threads
        num_Far_face          !- number of coupling face in this node
    integer,dimension(:),allocatable :: &
        size_Far_1D           !- size of independant component to couple
    real,dimension(:,:,:,:),allocatable :: &
        W_x1,W_x2,W_y1,W_y2,W_z1,W_z2
    real,dimension(:),allocatable :: &
        Vx_1D,Vy_1D,Vz_1D, &
        Txx_1D,Tyy_1D,Tzz_1D, &
        Tyz_1D,Txz_1D,Txy_1D

    integer ::  &
        siz_1D, &
        max_siz_1D, &
        indx_in_locl(2*SEIS_GEO,2*SEIS_GEO)    !- index in the local array
    integer,dimension(:),allocatable :: &
        fids                              !- id of the open files for coupling node
    integer,dimension(:,:,:),allocatable :: &
        indx_in_glob                      !- index in the global array of each node
    logical :: &
        flag_is,  & !- if coupling node
        flag_sid_node(2*SEIS_GEO) !- if coupled at this side
    logical,dimension(:,:),allocatable :: &
        flag_sid(:,:)
    character (len=SEIS_STRLEN) :: pnm_out

!-- loop var
    integer :: n,m,i,j,k,ii,nsid,n_i,n_j,n_k
    integer :: i1,i2,j1,j2,k1,k2,ntwin
    
!-------------------------------------------------------------------------------
! Program Entries:
!-------------------------------------------------------------------------------

    call get_conf_name(fnm_conf)
    call para_init(fnm_conf)
    call swmpi_init(fnm_conf)

!-------------------------------------------------------------------------------
!-- input parameters from standard input
!-------------------------------------------------------------------------------
    !write(*,*) "coupling index (x1 x2 y1 y2 z1 z2, 0 means no coupling):"
    write(*,*) "total far field simulation time steps:"
    read (*,*) ntwin_Far_field
    write(*,*) ntwin_Far_field

    write(*,*) "coupling at (x1)?:"
    read (*,*) gflag_x1
    write(*,*) gflag_x1
    write(*,*) "coupling index (x1):"
    read (*,*) gindx_x1
    write(*,*) gindx_x1
    write(*,*) "coupling filename (x1):"
    read (*,'(a)') gfnm_x1
    write(*,*) trim(gfnm_x1)

    write(*,*) "coupling at (x2)?:"
    read (*,*) gflag_x2
    write(*,*) gflag_x2
    write(*,*) "coupling index (x2):"
    read (*,*) gindx_x2
    write(*,*) gindx_x2
    write(*,*) "coupling filename (x2):"
    read (*,'(a)') gfnm_x2
    write(*,*) trim(gfnm_x2)

    write(*,*) "coupling at (y1)?:"
    read (*,*) gflag_y1
    write(*,*) gflag_y1
    write(*,*) "coupling index (y1):"
    read (*,*) gindx_y1
    write(*,*) gindx_y1
    write(*,*) "coupling filename (y1):"
    read (*,'(a)') gfnm_y1
    write(*,*) trim(gfnm_y1)

    write(*,*) "coupling at (y2)?:"
    read (*,*) gflag_y2
    write(*,*) gflag_y2
    write(*,*) "coupling index (y2):"
    read (*,*) gindx_y2
    write(*,*) gindx_y2
    write(*,*) "coupling filename (y2):"
    read (*,'(a)') gfnm_y2
    write(*,*) trim(gfnm_y2)

    write(*,*) "coupling at (z1)?:"
    read (*,*) gflag_z1
    write(*,*) gflag_z1
    write(*,*) "coupling index (z1):"
    read (*,*) gindx_z1
    write(*,*) gindx_z1
    write(*,*) "coupling filename (z1):"
    read (*,'(a)') gfnm_z1
    write(*,*) trim(gfnm_z1)

    write(*,*) "coupling at (z2)?:"
    read (*,*) gflag_z2
    write(*,*) gflag_z2
    write(*,*) "coupling index (z2):"
    read (*,*) gindx_z2
    write(*,*) gindx_z2
    write(*,*) "coupling filename (z2):"
    read (*,'(a)') gfnm_z2
    write(*,*) trim(gfnm_z2)

    write(*,*) "output dir:"
    read (*,'(a)') pnm_out
    write(*,*) trim(pnm_out)

!-- theoretical maximum coupling threads
    max_thread_couple=dims(1)*2+dims(2)*2+dims(3)*2

    allocate(fids         (max_thread_couple));
    allocate(size_Far_1D  (max_thread_couple));
    allocate(indx_in_glob (2*SEIS_GEO,2*SEIS_GEO,max_thread_couple))
    allocate(flag_sid     (           2*SEIS_GEO,max_thread_couple))

!-------------------------------------------------------------------------------
!-- create output files
!-------------------------------------------------------------------------------
    n=0   !- count of coupling node
    m=100
    max_siz_1D=0
    do n_i=0,dims(1)-1
    do n_j=0,dims(2)-1
    do n_k=0,dims(3)-1
        write(*,"(i10,2(i2),a,3(i2))") n_i,n_j,n_k, ' of ',dims
        call swmpi_change_fnm(n_i,n_j,n_k)
        call swmpi_set_gindx(n_i,n_j,n_k)

        m=m+1

        open(m+1,file=macdrp_get_fnm_coupling(pnm_out,n_i,n_j,n_k), &
                 status='unknown',form='unformatted')

        num_Far_face=0
        siz_1D=0
        flag_sid_node=.false.
        indx_in_locl(1:5:2,:)=1
        indx_in_locl(2:6:2,:)=0
        !-- for x,y,z overlap region
        indx_in_locl(2,1)=0; indx_in_locl(1,2)=ni+1
        indx_in_locl(4,3)=0; indx_in_locl(3,4)=nj+1
        !--
        if ( gflag_x1 .and. is_coupling_node(n_i,n_j,n_k,gindx_x1) ) then
           flag_sid_node(1)=.true.
           !-- index
           indx_in_locl(1,1)=gindx_x1(1,1)-n_i*ni
           indx_in_locl(2,1)=gindx_x1(2,1)-n_i*ni
           indx_in_locl(3,1)=max( n_j*nj+1,gindx_x1(1,2) )   - n_j*nj
           indx_in_locl(4,1)=min( (n_j+1)*nj,gindx_x1(2,2) ) - n_j*nj
           indx_in_locl(5,1)=max( n_k*nk+1,gindx_x1(1,3) )   - n_k*nk
           indx_in_locl(6,1)=min( (n_k+1)*nk,gindx_x1(2,3) ) - n_k*nk
           !-- size
           siz_1D  = siz_1D  &
                  + (indx_in_locl(2,1)-indx_in_locl(1,1)+1) &
                  * (indx_in_locl(4,1)-indx_in_locl(3,1)+1) &
                  * (indx_in_locl(6,1)-indx_in_locl(5,1)+1)
        end if
        if ( gflag_x2 .and. is_coupling_node(n_i,n_j,n_k,gindx_x2) ) then
           flag_sid_node(2)=.true.
           !-- index
           indx_in_locl(1,2)=gindx_x2(1,1)-n_i*ni
           indx_in_locl(2,2)=gindx_x2(2,1)-n_i*ni
           indx_in_locl(3,2)=max( n_j*nj+1,gindx_x2(1,2) )   - n_j*nj
           indx_in_locl(4,2)=min( (n_j+1)*nj,gindx_x2(2,2) ) - n_j*nj
           indx_in_locl(5,2)=max( n_k*nk+1,gindx_x2(1,3) )   - n_k*nk
           indx_in_locl(6,2)=min( (n_k+1)*nk,gindx_x2(2,3) ) - n_k*nk
           !-- size
           siz_1D  = siz_1D  &
                  + (indx_in_locl(2,2)-indx_in_locl(1,2)+1) &
                  * (indx_in_locl(4,2)-indx_in_locl(3,2)+1) &
                  * (indx_in_locl(6,2)-indx_in_locl(5,2)+1)
        end if
        if ( gflag_y1 .and. is_coupling_node(n_i,n_j,n_k,gindx_y1) ) then
           flag_sid_node(3)=.true.
           !-- index
           indx_in_locl(1,3)=max(max( n_i*ni+1,  gindx_y1(1,1) ) - n_i*ni, indx_in_locl(2,1)+1 ) !- remove points in x1
           indx_in_locl(2,3)=min(min( (n_i+1)*ni,gindx_y1(2,1) )-n_i*ni, indx_in_locl(1,2)-1 ) !- remove points in x2
           indx_in_locl(3,3)=gindx_y1(1,2)-n_j*nj
           indx_in_locl(4,3)=gindx_y1(2,2)-n_j*nj
           indx_in_locl(5,3)=max( n_k*nk+1,  gindx_y1(1,3) )   - n_k*nk
           indx_in_locl(6,3)=min( (n_k+1)*nk,gindx_y1(2,3) ) - n_k*nk
           !-- size
           siz_1D  = siz_1D  &
                  + (indx_in_locl(2,3)-indx_in_locl(1,3)+1) &
                  * (indx_in_locl(4,3)-indx_in_locl(3,3)+1) &
                  * (indx_in_locl(6,3)-indx_in_locl(5,3)+1)
        end if
        if ( gflag_y2 .and. is_coupling_node(n_i,n_j,n_k,gindx_y2) ) then
           flag_sid_node(4)=.true.
           !-- index
           indx_in_locl(1,4)=max(max( n_i*ni+1,  gindx_y2(1,1) ) - n_i*ni, indx_in_locl(2,1)+1 ) !- remove points in x1
           indx_in_locl(2,4)=min(min( (n_i+1)*ni,gindx_y2(2,1) )-n_i*ni, indx_in_locl(1,2)-1 ) !- remove points in x2
           indx_in_locl(3,4)=gindx_y2(1,2)-n_j*nj
           indx_in_locl(4,4)=gindx_y2(2,2)-n_j*nj
           indx_in_locl(5,4)=max( n_k*nk+1,  gindx_y2(1,3) )   - n_k*nk
           indx_in_locl(6,4)=min( (n_k+1)*nk,gindx_y2(2,3) ) - n_k*nk
           !-- size
           siz_1D  = siz_1D  &
                  + (indx_in_locl(2,4)-indx_in_locl(1,4)+1) &
                  * (indx_in_locl(4,4)-indx_in_locl(3,4)+1) &
                  * (indx_in_locl(6,4)-indx_in_locl(5,4)+1)
        end if
        if ( gflag_z1 .and. is_coupling_node(n_i,n_j,n_k,gindx_z1) ) then
           flag_sid_node(5)=.true.
           !-- index
           indx_in_locl(1,5)=max(max( n_i*ni+1,  gindx_z1(1,1) ) - n_i*ni, indx_in_locl(2,1)+1 ) !- remove points in x1
           indx_in_locl(2,5)=min(min( (n_i+1)*ni,gindx_z1(2,1) )-n_i*ni, indx_in_locl(1,2)-1 ) !- remove points in x2
           indx_in_locl(3,5)=max(max( n_j*nj+1,  gindx_z1(1,2) ) - n_j*nj, indx_in_locl(4,3)+1 ) !- remove points in y1
           indx_in_locl(4,5)=min(min( (n_j+1)*nj,gindx_z1(2,2) )-n_j*nj, indx_in_locl(3,4)-1 ) !- remove points in y2
           indx_in_locl(5,5)=gindx_z1(1,3)-n_k*nk
           indx_in_locl(6,5)=gindx_z1(2,3)-n_k*nk
           !-- size
           siz_1D  = siz_1D  &
                  + (indx_in_locl(2,5)-indx_in_locl(1,5)+1) &
                  * (indx_in_locl(4,5)-indx_in_locl(3,5)+1) &
                  * (indx_in_locl(6,5)-indx_in_locl(5,5)+1)
        end if
        if ( gflag_z2 .and. is_coupling_node(n_i,n_j,n_k,gindx_z2) ) then
           flag_sid_node(6)=.true.
           !-- index
           indx_in_locl(1,6)=max(max( n_i*ni+1,  gindx_z2(1,1) ) - n_i*ni, indx_in_locl(2,1)+1 ) !- remove points in x1
           indx_in_locl(2,6)=min(min( (n_i+1)*ni,gindx_z2(2,1) )-n_i*ni, indx_in_locl(1,2)-1 ) !- remove points in x2
           indx_in_locl(3,6)=max(max( n_j*nj+1,  gindx_z2(1,2) ) - n_j*nj, indx_in_locl(4,3)+1 ) !- remove points in y1
           indx_in_locl(4,6)=min(min( (n_j+1)*nj,gindx_z2(2,2) )-n_j*nj, indx_in_locl(3,4)-1 ) !- remove points in y2
           indx_in_locl(5,6)=gindx_z2(1,3)-n_k*nk
           indx_in_locl(6,6)=gindx_z2(2,3)-n_k*nk
           !-- size
           siz_1D  = siz_1D  &
                  + (indx_in_locl(2,6)-indx_in_locl(1,6)+1) &
                  * (indx_in_locl(4,6)-indx_in_locl(3,6)+1) &
                  * (indx_in_locl(6,6)-indx_in_locl(5,6)+1)
        end if

    !-- keep thread information if coupling node
        if (any(flag_sid_node)) then
           flag_is=.true.
           n=n+1
           fids(n)=m+1
           flag_sid(:,n)=flag_sid_node
           size_Far_1D(n)=siz_1D
           max_siz_1D = max(max_siz_1D,siz_1D)
           call conv_locl_glob(indx_in_locl,n_i,n_j,n_k,indx_in_glob(:,:,n))
        else
           flag_is=.false.
        end if

    !-- write index information
        write(m+1) flag_is,siz_1D,indx_in_locl,ntwin_Far_field/2

    !-- dimension size
        write(m+1) ni,nj,nk

     !-- close the file if not a coupling node
        if (.not. flag_is) close(m+1)
     
    end do !n_k
    end do !n_j
    end do !n_i

!-- reset max_thread_couple according to true value
    max_thread_couple = n

!-- allocate 1D var
    allocate( Vx_1D(max_siz_1D));
    allocate( Vy_1D(max_siz_1D));
    allocate( Vz_1D(max_siz_1D));
    allocate(Txx_1D(max_siz_1D));
    allocate(Tyy_1D(max_siz_1D));
    allocate(Tzz_1D(max_siz_1D));
    allocate(Tyz_1D(max_siz_1D));
    allocate(Txz_1D(max_siz_1D));
    allocate(Txy_1D(max_siz_1D));

!-------------------------------------------------------------------------------
!-- read in data and write out to node files
!-------------------------------------------------------------------------------
    if (gflag_x1) then
       allocate(W_x1(gindx_x1(1,1):gindx_x1(2,1),gindx_x1(1,2):gindx_x1(2,2),gindx_x1(1,3):gindx_x1(2,3),9))
       open(1,file=trim(gfnm_x1),status='old',form='unformatted')
       read(1) i1,i2,j1,j2,k1,k2,ntwin
       write(*,*) "in "//trim(gfnm_x1)//":"
       write(*,"(5i5,i9)") "    ",i1,i2,j1,j2,k1,k2,ntwin
       if (     i1/=gindx_x1(1,1) .or. i2/=gindx_x1(2,1) &
           .or. j1/=gindx_x1(1,2) .or. j2/=gindx_x1(2,2) &
           .or. k1/=gindx_x1(1,3) .or. k2/=gindx_x1(2,3) &
           .or. ntwin<ntwin_Far_field) then
           write(*,*) "header in "//trim(gfnm_x1)//" incoincident with input"
           stop 1
        end if
    end if
    if (gflag_x2) then
       allocate(W_x2(gindx_x2(1,1):gindx_x2(2,1),gindx_x2(1,2):gindx_x2(2,2),gindx_x2(1,3):gindx_x2(2,3),9))
       open(2,file=trim(gfnm_x2),status='old',form='unformatted')
       read(2) i1,i2,j1,j2,k1,k2,ntwin
       write(*,*) "in "//trim(gfnm_x2)//":"
       write(*,"(5i5,i9)") "    ",i1,i2,j1,j2,k1,k2,ntwin
       if (     i1/=gindx_x2(1,1) .or. i2/=gindx_x2(2,1) &
           .or. j1/=gindx_x2(1,2) .or. j2/=gindx_x2(2,2) &
           .or. k1/=gindx_x2(1,3) .or. k2/=gindx_x2(2,3) &
           .or. ntwin<ntwin_Far_field) then
           write(*,*) "header in "//trim(gfnm_x2)//" incoincident with input"
           stop 1
        end if
    end if
    if (gflag_y1) then
       allocate(W_y1(gindx_y1(1,1):gindx_y1(2,1),gindx_y1(1,2):gindx_y1(2,2),gindx_y1(1,3):gindx_y1(2,3),9))
       open(3,file=trim(gfnm_y1),status='old',form='unformatted')
       read(3) i1,i2,j1,j2,k1,k2,ntwin
       write(*,*) "in "//trim(gfnm_y1)//":"
       write(*,"(5i5,i9)") "    ",i1,i2,j1,j2,k1,k2,ntwin
       if (     i1/=gindx_y1(1,1) .or. i2/=gindx_y1(2,1) &
           .or. j1/=gindx_y1(1,2) .or. j2/=gindx_y1(2,2) &
           .or. k1/=gindx_y1(1,3) .or. k2/=gindx_y1(2,3) &
           .or. ntwin<ntwin_Far_field) then
           write(*,*) "header in "//trim(gfnm_y1)//" incoincident with input"
           stop 1
        end if
    end if
    if (gflag_y2) then
       allocate(W_y2(gindx_y2(1,1):gindx_y2(2,1),gindx_y2(1,2):gindx_y2(2,2),gindx_y2(1,3):gindx_y2(2,3),9))
       open(4,file=trim(gfnm_y2),status='old',form='unformatted')
       read(4) i1,i2,j1,j2,k1,k2,ntwin
       write(*,*) "in "//trim(gfnm_y2)//":"
       write(*,"(5i5,i9)") "    ",i1,i2,j1,j2,k1,k2,ntwin
       if (     i1/=gindx_y2(1,1) .or. i2/=gindx_y2(2,1) &
           .or. j1/=gindx_y2(1,2) .or. j2/=gindx_y2(2,2) &
           .or. k1/=gindx_y2(1,3) .or. k2/=gindx_y2(2,3) &
           .or. ntwin<ntwin_Far_field) then
           write(*,*) "header in "//trim(gfnm_y2)//" incoincident with input"
           stop 1
        end if
    end if
    if (gflag_z1) then
       allocate(W_z1(gindx_z1(1,1):gindx_z1(2,1),gindx_z1(1,2):gindx_z1(2,2),gindx_z1(1,3):gindx_z1(2,3),9))
       open(5,file=trim(gfnm_z1),status='old',form='unformatted')
       read(5) i1,i2,j1,j2,k1,k2,ntwin
       write(*,*) "in "//trim(gfnm_z1)//":"
       write(*,"(5i5,i9)") "    ",i1,i2,j1,j2,k1,k2,ntwin
       if (     i1/=gindx_z1(1,1) .or. i2/=gindx_z1(2,1) &
           .or. j1/=gindx_z1(1,2) .or. j2/=gindx_z1(2,2) &
           .or. k1/=gindx_z1(1,3) .or. k2/=gindx_z1(2,3) &
           .or. ntwin<ntwin_Far_field) then
           write(*,*) "header in "//trim(gfnm_z1)//" incoincident with input"
           stop 1
        end if
    end if
    if (gflag_z2) then
       allocate(W_z2(gindx_z2(1,1):gindx_z2(2,1),gindx_z2(1,2):gindx_z2(2,2),gindx_z2(1,3):gindx_z2(2,3),9))
       open(6,file=trim(gfnm_z2),status='old',form='unformatted')
       read(6) i1,i2,j1,j2,k1,k2,ntwin
       write(*,*) "in "//trim(gfnm_z2)//":"
       write(*,"(5i5,i9)") "    ",i1,i2,j1,j2,k1,k2,ntwin
       if (     i1/=gindx_z2(1,1) .or. i2/=gindx_z2(2,1) &
           .or. j1/=gindx_z2(1,2) .or. j2/=gindx_z2(2,2) &
           .or. k1/=gindx_z2(1,3) .or. k2/=gindx_z2(2,3) &
           .or. ntwin<ntwin_Far_field) then
           write(*,*) "header in "//trim(gfnm_z2)//" incoincident with input"
           stop 1
        end if
    end if

    do n=1,ntwin_Far_field
        write(*,*) "step:",n
    !-- read from global
        if (gflag_x1) read(1) W_x1
#ifdef DEBUG
        write(*,*) "after x1"
#endif
        if (gflag_x2) read(2) W_x2
#ifdef DEBUG
        write(*,*) "after x2"
#endif
        if (gflag_y1) read(3) W_y1
#ifdef DEBUG
        write(*,*) "after y1"
#endif
        if (gflag_y2) read(4) W_y2
#ifdef DEBUG
        write(*,*) "after y2"
#endif
        if (gflag_z1) read(5) W_z1
#ifdef DEBUG
        write(*,*) "after z1"
#endif
        if (gflag_z2) read(6) W_z2
#ifdef DEBUG
        write(*,*) "after z2"
#endif
    !-- write Far data
        do m=1,max_thread_couple
           ii=0
           do nsid=1,2*SEIS_GEO
              if (flag_sid(nsid,m)) then
                 do k=indx_in_glob(5,nsid,m),indx_in_glob(6,nsid,m)
                 do j=indx_in_glob(3,nsid,m),indx_in_glob(4,nsid,m)
                 do i=indx_in_glob(1,nsid,m),indx_in_glob(2,nsid,m)
                    ii=ii+1
                    if (nsid==1) then
                         Vx_1D(ii)=W_x1(i,j,k,1)
                         Vy_1D(ii)=W_x1(i,j,k,2)
                         Vz_1D(ii)=W_x1(i,j,k,3)
                        Txx_1D(ii)=W_x1(i,j,k,4)
                        Tyy_1D(ii)=W_x1(i,j,k,5)
                        Tzz_1D(ii)=W_x1(i,j,k,6)
                        Tyz_1D(ii)=W_x1(i,j,k,7)
                        Txz_1D(ii)=W_x1(i,j,k,8)
                        Txy_1D(ii)=W_x1(i,j,k,9)
                    end if
                    if (nsid==2) then
                         Vx_1D(ii)=W_x2(i,j,k,1)
                         Vy_1D(ii)=W_x2(i,j,k,2)
                         Vz_1D(ii)=W_x2(i,j,k,3)
                        Txx_1D(ii)=W_x2(i,j,k,4)
                        Tyy_1D(ii)=W_x2(i,j,k,5)
                        Tzz_1D(ii)=W_x2(i,j,k,6)
                        Tyz_1D(ii)=W_x2(i,j,k,7)
                        Txz_1D(ii)=W_x2(i,j,k,8)
                        Txy_1D(ii)=W_x2(i,j,k,9)
                    end if
                    if (nsid==3) then
                         Vx_1D(ii)=W_y1(i,j,k,1)
                         Vy_1D(ii)=W_y1(i,j,k,2)
                         Vz_1D(ii)=W_y1(i,j,k,3)
                        Txx_1D(ii)=W_y1(i,j,k,4)
                        Tyy_1D(ii)=W_y1(i,j,k,5)
                        Tzz_1D(ii)=W_y1(i,j,k,6)
                        Tyz_1D(ii)=W_y1(i,j,k,7)
                        Txz_1D(ii)=W_y1(i,j,k,8)
                        Txy_1D(ii)=W_y1(i,j,k,9)
                    end if
                    if (nsid==4) then
                         Vx_1D(ii)=W_y2(i,j,k,1)
                         Vy_1D(ii)=W_y2(i,j,k,2)
                         Vz_1D(ii)=W_y2(i,j,k,3)
                        Txx_1D(ii)=W_y2(i,j,k,4)
                        Tyy_1D(ii)=W_y2(i,j,k,5)
                        Tzz_1D(ii)=W_y2(i,j,k,6)
                        Tyz_1D(ii)=W_y2(i,j,k,7)
                        Txz_1D(ii)=W_y2(i,j,k,8)
                        Txy_1D(ii)=W_y2(i,j,k,9)
                    end if
                    if (nsid==5) then
                         Vx_1D(ii)=W_z1(i,j,k,1)
                         Vy_1D(ii)=W_z1(i,j,k,2)
                         Vz_1D(ii)=W_z1(i,j,k,3)
                        Txx_1D(ii)=W_z1(i,j,k,4)
                        Tyy_1D(ii)=W_z1(i,j,k,5)
                        Tzz_1D(ii)=W_z1(i,j,k,6)
                        Tyz_1D(ii)=W_z1(i,j,k,7)
                        Txz_1D(ii)=W_z1(i,j,k,8)
                        Txy_1D(ii)=W_z1(i,j,k,9)
                    end if
                    if (nsid==6) then
                         Vx_1D(ii)=W_z2(i,j,k,1)
                         Vy_1D(ii)=W_z2(i,j,k,2)
                         Vz_1D(ii)=W_z2(i,j,k,3)
                        Txx_1D(ii)=W_z2(i,j,k,4)
                        Tyy_1D(ii)=W_z2(i,j,k,5)
                        Tzz_1D(ii)=W_z2(i,j,k,6)
                        Tyz_1D(ii)=W_z2(i,j,k,7)
                        Txz_1D(ii)=W_z2(i,j,k,8)
                        Txy_1D(ii)=W_z2(i,j,k,9)
                    end if
                 end do
                 end do
                 end do
              end if
           end do
           write(fids(m)) (  Vx_1D(ii),ii=1,size_Far_1D(m) ), &
                          (  Vy_1D(ii),ii=1,size_Far_1D(m) ), &
                          (  Vz_1D(ii),ii=1,size_Far_1D(m) ), &
                          ( Txx_1D(ii),ii=1,size_Far_1D(m) ), &
                          ( Tyy_1D(ii),ii=1,size_Far_1D(m) ), &
                          ( Tzz_1D(ii),ii=1,size_Far_1D(m) ), &
                          ( Tyz_1D(ii),ii=1,size_Far_1D(m) ), &
                          ( Txz_1D(ii),ii=1,size_Far_1D(m) ), &
                          ( Txy_1D(ii),ii=1,size_Far_1D(m) )
        end do !- m
    end do !- n

!-- close files
    do m=1,max_thread_couple
       close(fids(m))
    end do
    if (gflag_x1) close(1) 
    if (gflag_x2) close(2)
    if (gflag_y1) close(4)
    if (gflag_y2) close(5)
    if (gflag_z1) close(6)
    if (gflag_z2) close(7)

!-- deallocate
    if (allocated( Vx_1D)) deallocate( Vx_1D)
    if (allocated( Vy_1D)) deallocate( Vy_1D)
    if (allocated( Vz_1D)) deallocate( Vz_1D)
    if (allocated(Txx_1D)) deallocate(Txx_1D)
    if (allocated(Tyy_1D)) deallocate(Tyy_1D)
    if (allocated(Tzz_1D)) deallocate(Tzz_1D)
    if (allocated(Tyz_1D)) deallocate(Tyz_1D)
    if (allocated(Txz_1D)) deallocate(Txz_1D)
    if (allocated(Txy_1D)) deallocate(Txy_1D)
    if (allocated(W_x1)) deallocate(W_x1)
    if (allocated(W_x2)) deallocate(W_x2)
    if (allocated(W_y1)) deallocate(W_y1)
    if (allocated(W_y2)) deallocate(W_y2)
    if (allocated(W_z1)) deallocate(W_z1)
    if (allocated(W_z2)) deallocate(W_z2)

!===============================================================================
contains
!===============================================================================

!-------------------------------------------------------------------------------
function is_coupling_node(n_i,n_j,n_k,indx) result(is_not)
!-------------------------------------------------------------------------------
    integer :: n_i,n_j,n_k
    integer,dimension(2,SEIS_GEO) :: indx
    logical :: is_not
    is_not=.false.
    if (      indx(1,1)<= (n_i+1)*ni .and. indx(2,1)>= n_i*ni+1 &
        .and. indx(1,2)<= (n_j+1)*nj .and. indx(2,2)>= n_j*nj+1 &
        .and. indx(1,3)<= (n_k+1)*nk .and. indx(2,3)>= n_k*nk+1 ) then
       is_not=.true.
    end if
!-------------------------------------------------------------------------------
end function is_coupling_node
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine conv_locl_glob(indx_in_locl,n_i,n_j,n_k,indx_in_glob)
!-------------------------------------------------------------------------------
    integer :: n_i,n_j,n_k
    integer ::                               &
        indx_in_locl(2*SEIS_GEO,2*SEIS_GEO), & !- index in the local array
        indx_in_glob(2*SEIS_GEO,2*SEIS_GEO)    !- index in the global array

    indx_in_glob(1:2,:)=indx_in_locl(1:2,:)+n_i*ni
    indx_in_glob(3:4,:)=indx_in_locl(3:4,:)+n_j*nj
    indx_in_glob(5:6,:)=indx_in_locl(5:6,:)+n_k*nk

!-------------------------------------------------------------------------------
end subroutine conv_locl_glob
!-------------------------------------------------------------------------------

end program tool_coupling

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
