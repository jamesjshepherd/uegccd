
      Subroutine DiagR(Matrix,Evals,Evecs,N)
      Use Precision
      Implicit None
      Integer,        Intent(In)  :: N
      Real (Kind=pr), Intent(In)  :: Matrix(N,N)
      Real (Kind=pr), Intent(Out) :: Evecs(N,N), Evals(N)
      Real (Kind=pr), Allocatable :: Work(:)
      Integer :: NBlocks, LWork, Info, ILAENV

!==========================================!
!  Non-destructive wrapper to DSyEv.       !
!  Eigenvalues sorted in ascending order.  !
!==========================================!

! Find LWork and Allocate
      NBlocks = ILAENV(1,"DSYTRD","U",N,-1,-1,-1) + 2
      LWork   = N*NBlocks
      Allocate(Work(LWork),  Stat=Info)
      If(Info /= 0) Stop "Could not allocate in DiagR"

! Diagonalize!
      Evecs = Matrix
      Call DSyEv("V","U",N,Evecs,N,Evals,Work,LWork,Info)
      If(Info /= 0) Stop "Error in LAPack!"

! Deallocate and exit safely
      DeAllocate(Work,  Stat=Info)
      If(Info /= 0) Stop "Could not deallocate in DiagR"

      Return
      End Subroutine DiagR

