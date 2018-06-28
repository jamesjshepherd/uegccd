
      Function AllCase(String)

      Implicit None
      Character (len=*) :: String
      Character (len=*) :: AllCase
      Integer   :: I, J

!==============================================================!
!  Given a string, returns it with all letters in upper case.  !
!  Here, we've assumed that the lower case letters in the      !
!  ASCII sequence run from 97 to 122 and the upper case run    !
!  from 65 to 90.  We could put in more logic to handle other  !
!  possibilities, but let's not worry about that right now.    !
!==============================================================!

      AllCase = String
      Do I = 1,Len(String)
       J = IAChar(String(I:I))
       If(97 <= J .and. 122 >= J) AllCase(I:I) = AChar(J-32)
      End Do

      Return
      End Function AllCase

