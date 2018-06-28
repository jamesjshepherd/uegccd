
   Module Precision

   Implicit None
   Integer, Parameter :: dp = Selected_Real_Kind(15,307)
   Integer, Parameter :: pr = dp

!===================================================!
!  The precision we work in corresponds to double   !
!  precision, according to Paul; he gets this from  !
!  http://www.nsc.liu.se/~boein/f77to90/c13.html    !
!                                                   !
!  We've chosen the kind that gives us at least 15  !
!  decimal places and exponents in [-307,307].      !
!===================================================!

   End Module Precision


