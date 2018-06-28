
   Module DIISCC

   Use Precision
   Use Constants
   Implicit None
   Real (Kind=pr), Allocatable :: OldT2abab(:,:,:,:), R2abab(:,:,:,:)
   Real (Kind=pr), Allocatable :: OldT2abba(:,:,:,:), R2abba(:,:,:,:)

   Contains

      Subroutine SetUpDIISCC(O,V,N,NDIIS)
      Implicit None
      Integer, Intent(In)  :: O, V, N
      Integer, Intent(Out) :: NDIIS
      Integer, Parameter   :: Max_N_DIIS = 3
      Integer :: IAlloc

      Do NDIIS = Max_N_DIIS,1,-1
      Allocate(OldT2abab(1:O,1:O,V:N,NDIIS),                &
               R2abab(1:O,1:O,V:N,NDIIS),                   &
               OldT2abba(1:O,1:O,V:N,NDIIS),                &
               R2abba(1:O,1:O,V:N,NDIIS),                   &
               Stat=IAlloc)
        If(IAlloc == 0) Exit
      End Do

      If( .not. Allocated(OldT2abab)) Stop "Could not allocate for DIIS"

      Return
      End Subroutine SetUpDIISCC






      Subroutine GetTFromDIIS(T2,OldT2,R2,NOcc,NBF,NDIIS,NIter)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF, NDIIS, NIter
      Real (Kind=pr), Intent(Out)   :: T2(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: OldT2(NOcc,NOcc,NOcc+1:NBF,NDIIS)
      Real (Kind=pr), Intent(InOut) :: R2(NOcc,NOcc,NOcc+1:NBF,NDIIS)
! Local variables
      Real (Kind=pr) :: C(NDIIS), V(NDIIS+1), M(NDIIS+1,NDIIS+1)
      Integer :: IPiv(2*NDIIS+2)
      Integer :: I, J, A, B, P, Q, Info
      Integer, Parameter :: IWait = 5

!===============================================================!
!  Forms the DIIS guess for T, which we insert to calculate G.  !
!  Remember that DIIS writes                                    !
!     T* ~ c(i) T(i)                                            !
!  where the sum goes over the m previous estimates for T.  We  !
!  have the constraint                                          !
!     Sum(c(i)) = 1                                             !
!  and we get the c(i) by minimizing a residual, effectively:   !
!     r(i) = H.T(i) - G[T(i)]                                   !
!     R = c(i) r(i)                                             !
!  then minimize                                                !
!     R*R = Sum(c(i)* c(j) B(i,j))                              !
!     B(i,j) = r(i)*.r(j)                                       !
!                                                               !
!  Try the real case:                                           !
!     L = Sum[c(i) c(j) B(ij)] + 2 y (Sum[c(i)] - 1)            !
!  Then                                                         !
!     dL/dy = 0    => Sum(c(i)) = 1                             !
!     dL/dc(i) = 0 => 2 Sum(B(ij) cj) + 2 y = 0.                !
!  We can organize this as a matrix equation:                   !
!     [B(i,j)   1(j)] [c(j)] = [0(j)]                           !
!     [1(i)       0 ] [ y  ] = [ 1  ]                           !
!  which we define as M.C = V.                                  !
!===============================================================!

! Zero out C.  Set C(NDIIS) = 1 so that we can do this right if DIIS is turned off
      C = Zero
      C(NDIIS) = One

! Build the DIIS matrix M and vector V
      If(NIter >= NDIIS + IWait) Then
        V = Zero
        M = Zero
! Build the B-matrix part
        Do P = 1,NDIIS
        Do Q = 1,NDIIS
          Do I = 1,NOcc
          Do J = 1,NOcc
          Do A = NOcc+1,NBF
            M(P,Q) = M(P,Q) + R2(I,J,A,P)*R2(I,J,A,Q)
          End Do
          End Do
          End Do
        End Do
        End Do
! Build the one-vector part
        M(1:NDIIS,1:) = One
        M(1,1::NDIIS) = One
        V(NDIIS+1)    = One
! Solve for C
        Call DGeSV(NDIIS+1,1,M,NDIIS+1,IPiv,V,NDIIS+1,Info)
        If(Info == 0) C = V(1:NDIIS)

! Get the new T2
        T2 = Zero
        Do P = 1,NDIIS
          T2 = T2 + C(P)*OldT2(:,:,:,P)
        End Do
      Else
! We're not doing DIIS in this case
        T2 = OldT2(:,:,:,NDIIS)
      End If

! Update the residuals.  We'll need to build R2(NDIIS)
      R2(:,:,:,1:NDIIS-1) = R2(:,:,:,2:NDIIS)

      Return
      End Subroutine GetTFromDIIS






      Subroutine ShutDownDIISCC
      Implicit None
      Integer IAlloc
      DeAllocate(OldT2abab, OldT2abba, R2abab, R2abba, Stat=IAlloc)
      If(IAlloc /= 0) Stop "Could not deallocate for DIIS"
      Return
      End Subroutine ShutDownDIISCC

   End Module DIISCC

