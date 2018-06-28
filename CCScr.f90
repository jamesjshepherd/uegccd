
   Module CCScr

   Use Precision
   Implicit None
   Real (Kind=pr), Allocatable :: G2aaaa(:,:,:)
   Real (Kind=pr), Allocatable :: G2abab(:,:,:)
   Real (Kind=pr), Allocatable :: G2abba(:,:,:)
   Real (Kind=pr), Allocatable :: Joo(:), Jvv(:)

   Contains

      Subroutine SetUpCCScr(O,V,N)
      Implicit None
      Integer, Intent(In) :: O, V, N
      Integer :: IAlloc
      Allocate(G2aaaa(1:O,1:O,V:N),     &
               G2abab(1:O,1:O,V:N),     &
               G2abba(1:O,1:O,V:N),     &
               Joo(1:O), Jvv(V:N),      &
               Stat=IAlloc)
      If(IAlloc /= 0) Stop 'Could not allocate in CCScr'
      Return
      End Subroutine SetUpCCScr






      Subroutine ShutDownCCScr
      Implicit None
      Integer :: IAlloc
      DeAllocate(G2aaaa, G2abab, G2abba, Joo, Jvv, Stat=IAlloc)
      If(IAlloc /= 0) Stop 'Could not deallocate in CCScr'
      Return
      End Subroutine ShutDownCCScr

   End Module CCScr

