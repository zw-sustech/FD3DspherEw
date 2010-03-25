program test

    implicit none
    interface
        function bessi0_s(x) result(y)
            real :: x
            real :: y
        end function bessi0_s
    end interface

    integer,parameter :: NSRCEXT=4 !- for Kaiser and sinc function
    !real,parameter    :: KAISERB=4.14
    real,parameter    :: KAISERB=6.31
    real,parameter :: PI=3.1415926, SEIS_ZERO=1.0e-25
    real :: b0,x0,y0,z0
    integer :: Li
    real,dimension(-NSRCEXT:NSRCEXT) :: DNx, KAx
                
    
    b0=bessi0_s(KAISERB)

    !-- source position shift
        !x0=0.5
        read(*,*) x0
    !-- discrete delta function value
        do Li=-NSRCEXT,NSRCEXT
            if ( abs(Li-x0)<SEIS_ZERO ) then
                DNx(Li)=1.0
            else
                DNx(Li)=sin(PI*(Li-x0))/(PI*(Li-x0))
            end if
            !if ( Li+x0 <= NSRCEXT ) then
            if ( 1.0-(Li-x0)**2/NSRCEXT**2 > 0.0 ) then
                KAx(Li)=bessi0_s( KAISERB * sqrt( 1.0-(Li-x0)**2/NSRCEXT**2) ) / b0
            else
                KAx(Li)=0.0
            end if
            write(10,*) DNx(Li)
            write(11,*) KAx(Li)
            write(12,*) DNx(Li)*KAx(Li)
        end do

end program test
