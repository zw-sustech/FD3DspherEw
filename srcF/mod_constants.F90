!******************************************************************************!
!*  This module declares parameters used in other modules                     *!
!*                                                                            *!
!*  Author: Wei ZHANG     Email: zhangw.pku@gmail.com                         *!
!*  Copyright (C) Wei ZHANG, 2006. All Rights Reserved.                       *!
!******************************************************************************!

! $Date$
! $Revision$
! $LastChangedBy$

!-----------------------------------------------------------------------------
module constants_mod
!-----------------------------------------------------------------------------

implicit none

INTEGER, PARAMETER :: DP = KIND(1.0D0)

#ifdef DataTypeDouble
INTEGER, PARAMETER :: SP = KIND(1.0D0)
#else
INTEGER, PARAMETER :: SP = KIND(1.0)
#endif

real(DP),parameter ::        &
     PI=3.141592653589793_DP
real(SP),parameter ::        &
     SEIS_ZERO=1.0e-20_SP,   &
     SEIS_EQUAL=1.0e-5_SP,   &
     SEIS_NAN=1.0e+15_SP,    &
     SEIS_INF=1.0e+25_SP

integer,parameter ::    &
     SEIS_GEO=3,        &
     SEIS_STRLEN=132

end module constants_mod

! vim:ft=fortran:ts=4:sw=4:nu:et:ai:
