
   Module CCD

   Use Precision
   Use Constants
   Use CCScr    ! Holds G2abab, G2abba, Joo, Jvv
   Use DIISCC   ! Holds OldT2abab, OldT2abba, R2abab, R2abba
   Use HEG
   Implicit None
   Real (Kind=pr), Parameter :: DenomFactor = 2.0_pr

   Contains

      Subroutine DrvCCD(Fock,T2aaaa,T2abab,T2abba,NOcc,NBF,ESCF,ECorr,     &
                        DoRings,DoXRings,DoLadders,DoMosaics,              &
                        IRangeRing,IRangeXRing,IRangeLadder,IRangeMosaic,  &
                        IRangeDriverDirect,IRangeDriverExchange,IRangeEnergy,                       &
                        IRangeLinRings,IRangeQuadRings,IRangeDirectRings,IRangeExchangeRings, &
                        IRangeLinLadders,IRangeQuadLadders,IRangeDirectLadders,IRangeExchangeLadders)
      Use HEG, only: FindTol
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF
      Real (Kind=pr), Intent(In)    :: ESCF
      Real (Kind=pr), Intent(In)    :: Fock(NBF)
      Real (Kind=pr), Intent(Out)   :: ECorr
      Logical,        Intent(In)    :: DoRings, DoXRings, DoLadders, DoMosaics
      Integer,        Intent(In)    :: IRangeRing, IRangeXRing, IRangeLadder, IRangeMosaic
      Integer,        Intent(In)    :: IRangeDriverDirect, IRangeDriverExchange, IRangeEnergy
      Integer,        Intent(In)    :: IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings  
      Integer,        Intent(In)    :: IRangeLinLadders, IRangeQuadLadders, IRangeDirectLadders, IRangeExchangeLadders
      Real (Kind=pr), Intent(InOut) :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr)                :: TolMax = 1.0E-8_pr
      Real (Kind=pr)                :: FailRatio = 100
      Integer,        Parameter     :: MaxIter = 1000
      Logical :: Fail
      Integer :: NIter, NDIIS
      Real (Kind=pr) :: dT, MaxRes

!========================================================!
!  This routine does the entire CCD calculation for us.  !
!  We take in the MBPT(2) results as input and write     !
!  them out as the first CCD iteration even though we    !
!  don't actually repeat them.                           !
!                                                        !
!  Certain terms of the CCD equations can be omitted in  !
!  what I think is an obvious way.                       !
!                                                        !
!  We use DIIS to converge this better.                  !
!========================================================!

      Call FindTol(TolMax,FailRatio)

      Write(6,*) 'Doing CCD...'
! Allocate space for everything.
      Call SetUpCCScr(NOcc,NOcc+1,NBF)
      Call SetUpDIISCC(NOcc,NOcc+1,NBF,NDIIS)


! Initialize variables and write the MP2 energy out again
      OldT2abab = Zero; OldT2abab(:,:,:,NDIIS) = T2abab; R2abab = Zero
      OldT2abba = Zero; OldT2abba(:,:,:,NDIIS) = T2abba; R2abba = Zero
      NIter = 1
      dT = Max(MaxVal(Abs(T2abab)),MaxVal(Abs(T2abba)))
      Call CCEnergy(T2aaaa,T2abab,T2abba,ECorr,NOcc,NBF,IRangeEnergy)
      Open(7,File='Output',Position='Append')
      Write(7,1010)
      Write(7,1020)
      If(DoRings)   Write(7,1100)
      If(DoXRings)  Write(7,1110)
      If(DoLadders) Write(7,1120)
      If(DoMosaics) Write(7,1130)
      Write(7,1020)
      Write(7,1030)
      Write(7,1040) ECorr,NIter,dT


! Start iterating!
      Do While(dT >= TolMax .and. NIter < MaxIter)
! Given OldT2 and R2, update T2 and R2
        Call GetTFromDIIS(T2abab,OldT2abab,R2abab,NOcc,NBF,NDIIS,NIter)
        Call GetTFromDIIS(T2abba,OldT2abba,R2abba,NOcc,NBF,NDIIS,NIter)
        T2aaaa = T2abab + T2abba
! Given T2 from DIIS, calculate G2.
        Call GetG2Drivers(G2abab,G2abba,NOcc,NBF,IRangeDriverDirect,IRangeDriverExchange)
        If(DoMosaics) Call GetG2Mosaics(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeMosaic)
        If(DoLadders) Call GetG2Ladders(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeLadder, &
                                         IRangeLinLadders, IRangeQuadLadders, IRangeDirectLadders, IRangeExchangeLadders)
        If(DoRings)   Call GetG2Rings(  G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeRing,  &
                                        IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings)
        If(DoXRings)  Call GetG2XRings( G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeXRing, &
                                        IRangeLinRings,IRangeQuadRings,IRangeDirectRings,IRangeExchangeRings)
        
! Calculate the residual G[T]-HT
        Call GetRes(T2abab,G2abab,R2abab,Fock,NOcc,NBF,NDIIS)
        Call GetRes(T2abba,G2abba,R2abba,Fock,NOcc,NBF,NDIIS)
! Solve the CC equations, update OldT2, and get the energy
        Call SolveCC(T2abab,G2abab,Fock,NOcc,NBF)
        Call SolveCC(T2abba,G2abba,Fock,NOcc,NBF)
        T2aaaa = T2abab + T2abba
        Call CnvrgCC(T2abab,OldT2abab,T2abba,OldT2abba,dT,NOcc,NBF,NIter,NDIIS)
        Call CCEnergy(T2aaaa,T2abab,T2abba,ECorr,NOcc,NBF,IRangeEnergy)
        Write(7,1040) ECorr, NIter, dT
        Write(6,1040) ECorr, NIter, dT
      End Do


! Check the resisdual
      Call GetG2Drivers(G2abab,G2abba,NOcc,NBF,IRangeDriverDirect,IRangeDriverExchange)
      If(DoMosaics) Call GetG2Mosaics(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeMosaic)
      If(DoLadders) Call GetG2Ladders(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeLadder, &
                                         IRangeLinLadders, IRangeQuadLadders, IRangeDirectLadders, IRangeExchangeLadders)
      If(DoRings)   Call GetG2Rings(  G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeRing, &
                                        IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings)
      If(DoXRings)  Call GetG2XRings( G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRangeXRing, &
                                        IRangeLinRings,IRangeQuadRings,IRangeDirectRings,IRangeExchangeRings)
      Call GetRes(T2abab,G2abab,R2abab,Fock,NOcc,NBF,NDIIS)
      Call GetRes(T2abba,G2abba,R2abba,Fock,NOcc,NBF,NDIIS)
      MaxRes =  Max(MaxVal(Abs(R2abab(:,:,:,NDIIS))),    &
                    MaxVal(Abs(R2abba(:,:,:,NDIIS))))
      Fail = MaxRes > FailRatio*TolMax


! Finish the ouptut and shut down
      Write(7,1020)
      If(Fail) Then
        Write(7,2000)
      Else
        Write(7,1050) NIter
        Write(7,1060) ECorr
        Write(7,1070) ECorr+ESCF
      End If
      Write(7,1090) MaxRes
      Write(7,1000)
      Close(7)
      Call ShutDownCCScr
      Call ShutDownDIISCC


1000  Format(14x,'**************************************************')
1010  Format(14X,'*               CCD summary follows              *')
1020  Format(14x,'*------------------------------------------------*')
1030  Format(14X,'*  CCD Energy      Iteration    Biggest change   *')
1040  Format(14X,'* ',F15.10,2x,I5,7X,F14.10,4x,'*')
1050  Format(14x,'* CCD has converged in ',I3,' iterations',12x,'*')
1060  Format(14x,'* Ec(CCD) is ',10x,F18.12,' a.u.   *')
1070  Format(14x,'* Final CCD Energy is  ',F18.12,' a.u.   *')
1090  Format(14x,'* Max CCD Residual is  ',4x,ES18.12,'    *')
1100  Format(14x,'* Including ring diagrams...                     *')
1110  Format(14x,'* Including crossed ring diagrams...             *')
1120  Format(14x,'* Including ladder diagrams...                   *')
1130  Format(14x,'* Including mosaic diagrams...                   *')
2000  Format(14x,'* Final CCD Energy is   DID NOT CONVERGE         *')

      Return
      End Subroutine DrvCCD






      Subroutine CCEnergy(T2aaaa,T2abab,T2abba,ECorr,NOcc,NBF,IRange)
      Implicit None
      Integer,        Intent(In)  :: NOcc,NBF
      Real (Kind=pr), Intent(In)  :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)  :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)  :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(Out) :: ECorr
      Real (Kind=pr) :: V_ijab, V_ijba, Eaaaa, Eabab, Eabba, dE
      Integer        :: I, J, A, B
      Integer,        Intent(In)    :: IRange

!===========================================!
!  This routine calculates the CCD energy.  !
!===========================================!

! THis needs to be fixed for momentum symmetry
      ECorr = Zero
      Eaaaa = Zero
      Eabab = Zero
      Eabba = Zero
      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle
        V_ijab = ERI(A,B,I,J,IRange)
        V_ijba = ERI(A,B,J,I,IRange)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
!       Eaaaa = Eaaaa + T2aaaa(I,J,A)*(V_ijab - V_ijba)
!       Eabab = Eabab + T2abab(I,J,A)*V_ijab
!       Eabba = Eabba - T2abba(I,J,A)*V_ijba
        dE = V_ijab*(T2aaaa(I,J,A) + T2abab(I,J,A))    &
           - V_ijba*(T2aaaa(I,J,A) + T2abba(I,J,A))
        ECorr = ECorr + F12*dE
      End Do
      End Do
      End Do
!     Eaaaa = Eaaaa*F12
!     Eabab = Eabab*F12
!     Eabba = Eabba*F12
!     ECorr = ECorr + Eaaaa + Eabab + Eabba

      Return
      End Subroutine CCEnergy






      Subroutine SolveCC(T2,G2,Fock,NOcc,NBF)
      Implicit None
      Integer,        Intent(In)  :: NOcc, NBF
      Real (Kind=pr), Intent(In)  :: Fock(NBF)
      Real (Kind=pr), Intent(In)  :: G2(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(Out) :: T2(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr) :: Denom
      Integer        :: I, J, A, B

!=============================!
!  Solves H.T = G.  Trivial.  !
!=============================!

      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle
        Denom = Fock(I) + Fock(J) - Fock(A) - Fock(B)
!       T2(I,J,A) = G2(I,J,A)/Denom
        T2(I,J,A) = (G2(I,J,A)/Denom - (One-DenomFactor)*T2(I,J,A))/DenomFactor
      End Do
      End Do
      End Do

      Return
      End Subroutine SolveCC






      Subroutine GetRes(T2,G2,R2,Fock,NOcc,NBF,NDIIS)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF, NDIIS
      Real (Kind=pr), Intent(In)    :: Fock(NBF)
      Real (Kind=pr), Intent(In)    :: T2(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: G2(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: R2(NOcc,NOcc,NOcc+1:NBF,NDIIS)
      Real (Kind=pr) :: Rijab
      Integer :: I, J, A, B

!===============================!
!  Computes the residual HT-G.  !
!===============================!

      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle
        Rijab = -G2(I,J,A) + (Fock(I) + Fock(J) - Fock(A) - Fock(B))*T2(I,J,A)
        R2(I,J,A,NDIIS) = Rijab
      End Do
      End Do
      End Do

      Return
      End Subroutine GetRes






      Subroutine CnvrgCC(T2abab,OldT2abab,T2abba,OldT2abba,dT,NOcc,NBF,NIter,NDIIS)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF, NDIIS
      Integer,        Intent(InOut) :: NIter
      Real (Kind=pr), Intent(Out)   :: dT
      Real (Kind=pr), Intent(In)    :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: OldT2abab(NOcc,NOcc,NOcc+1:NBF,NDIIS)
      Real (Kind=pr), Intent(InOut) :: OldT2abba(NOcc,NOcc,NOcc+1:NBF,NDIIS)
      Integer, Parameter :: MaxIter = 10000
      Real (Kind=pr) :: dTabab, dTabba

!======================================!
!  Update the DIIS amplitudes OldT2*.  !
!======================================!

! Update the iteration count
      NIter = NIter + 1

! Move the DIIS amplitudes
      OldT2abab(:,:,:,1:NDIIS-1) = OldT2abab(:,:,:,2:NDIIS)
      OldT2abba(:,:,:,1:NDIIS-1) = OldT2abba(:,:,:,2:NDIIS)

! Get the change in T
      dTabab = MaxVal(Abs(T2abab - OldT2abab(:,:,:,NDIIS)))
      dTabba = MaxVal(Abs(T2abba - OldT2abba(:,:,:,NDIIS)))
      dT = Max(dTabab, dTabba)

! Update the most recent DIIS amplitudes
      OldT2abab(:,:,:,NDIIS) = T2abab
      OldT2abba(:,:,:,NDIIS) = T2abba
      If(NIter == MaxIter) Stop 'CCD equations did not converge'

      Return
      End Subroutine CnvrgCC






      Subroutine GetG2Drivers(G2abab,G2abba,NOcc,NBF,IRangeDirect,IRangeExchange)
      Implicit None
      Integer,        Intent(In)  :: NOcc, NBF, IRangeDirect, IRangeExchange
      Real (Kind=pr), Intent(Out) :: G2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(Out) :: G2abba(NOcc,NOcc,NOcc+1:NBF)
! Local variables
      Integer :: I, J, A, B
      Real (Kind=pr) :: V_ijab, V_ijba

!==================================!
!  Build the driving terms first.  !
!==================================!

      G2abab = Zero
      G2abba = Zero
      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle
        V_ijab = ERI(I,J,A,B,IRangeDirect)
        V_ijba = ERI(I,J,B,A,IRangeExchange)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        G2abab(I,J,A) =  V_ijab
        G2abba(I,J,A) = -V_ijba
      End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Drivers






      Subroutine GetG2Mosaics(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRange)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF
      Real (Kind=pr), Intent(In)    :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abba(NOcc,NOcc,NOcc+1:NBF)
      Integer,        Intent(In)    :: IRange
! Local variables
      Integer :: I, J, A, B, K, L, C, D
      Real (Kind=pr) :: V_cdkl, V_cdlk

!===========================================================!
!  Handle the mosaic terms.  These are                      !
!   dG_ij^ab = -1/2 vbar^{kl}_{cd} t_{kl}^{ad} t_{ij}^{cb}  !
!              -1/2 vbar^{kl}_{cd} t_{kl}^{cb} t_{ij}^{ad}  !
!              -1/2 vbar^{kl}_{cd} t_{il}^{cd} t_{kj}^{ab}  !
!              -1/2 vbar^{kl}_{cd} t_{kj}^{cd} t_{il}^{ab}  !
!  They're particularly simple - we only have to build the  !
!  diagonal intermediates in this case.                     !
!                                                           !
!  Let's see that explicitly.  Consider first Jvv.          !
!  We would have                                            !
!     J^a_c = t_kl^ad v^kl_cd                               !
!  and then build                                           !
!     dG_ij^ab = t_ij^cb J^a_c                              !
!  But we know that i+j+a+b = 0, and we need i+j+c+b = 0,   !
!  which means that a=c => we need only the diagonal.       !
!-----------------------------------------------------------!
!  Start by building intermediates.                         !
!  It sure looks like the Joo and Jvv are identical!        !
!===========================================================!

      Joo = Zero
      Jvv = Zero
      Do K = 1,NOcc
      Do L = 1,NOcc
      Do C = NOcc+1,NBF
        D = FindIndex(K,L,C)
        If(D <= NOcc) Cycle
        V_cdkl = ERI(C,D,K,L,dummy_flag=IRange)
        V_cdlk = ERI(C,D,L,K,dummy_flag=IRange)
        If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        Jvv(C) = Jvv(C) + V_cdkl*(T2aaaa(K,L,C) + T2abab(K,L,C))               &
                        - V_cdlk*(T2aaaa(K,L,C) + T2abba(K,L,C))
        Joo(K) = Joo(K) + V_cdkl*(T2aaaa(K,L,C) + T2abab(K,L,C))               &
                        - V_cdlk*(T2aaaa(K,L,C) + T2abba(K,L,C))
      End Do
      End Do
      End Do

      
! The intermediates are done, so contract with T2
      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle
        G2abab(I,J,A) = G2abab(I,J,A) - F12*(Jvv(A)*T2abab(I,J,A) + Jvv(B)*T2abab(I,J,A))
        G2abba(I,J,A) = G2abba(I,J,A) - F12*(Jvv(A)*T2abba(I,J,A) + Jvv(B)*T2abba(I,J,A))
        G2abab(I,J,A) = G2abab(I,J,A) - F12*(Joo(I)*T2abab(I,J,A) + Joo(J)*T2abab(I,J,A))
        G2abba(I,J,A) = G2abba(I,J,A) - F12*(Joo(I)*T2abba(I,J,A) + Joo(J)*T2abba(I,J,A))
      End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Mosaics







      Subroutine GetG2Ladders(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRange, &
                              IRangeLinLadders, IRangeQuadLadders, IRangeDirectLadders, IRangeExchangeLadders)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF
      Real (Kind=pr), Intent(In)    :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abba(NOcc,NOcc,NOcc+1:NBF)
      Integer,        Intent(In)    :: IRange
      Integer,        Intent(In)    :: IRangeLinLadders, IRangeQuadLadders, IRangeDirectLadders, IRangeExchangeLadders
! Local variables
      Integer :: I, J, A, B, K, L, C, D
      Real (Kind=pr) :: J2_ijkl, J3_ijkl, V_cdkl, V_cdlk, V_cdab, V_cdba

!=================================================================!
!  Gets the ladder contributions here.  Specialized for the HEG.  !
!  As vectors, we have A + B = I + J = K + L = C + D.             !
!=================================================================!

      Do I = 1,NOcc
      Do J = 1,NOcc

! Build the intermediates for each IJKL => hole-hole ladder, ladder T.V.T
! We'll have J2_ijkl which goes with T2abab(K,L,A,B) in G2abab(I,J,A,B) and with T2abba(K,L,A,B) in G2abba(I,J,A,B)
! We'll have J3_ijkl which goes with T2abba(K,L,A,B) in G2abab(I,J,A,B) and with T2abab(K,L,A,B) in G2abba(I,J,A,B)
        Do K = 1,NOcc
          L = FindIndex(I,J,K)
          If(L > NOcc .or. L <= 0) Cycle
          J2_ijkl =  ERI(I,J,K,L,max(IRange,IRangeDirectLadders,IRangeLinLadders))              ! For T2abab in G2abab and T2abba in G2abba
          J3_ijkl = -ERI(I,J,L,K,max(IRange,IRangeExchangeLadders,IRangeLinLadders))              ! For T2abba in G2abab and T2abab in G2abba
          If(Min(J2_ijkl,-J3_ijkl) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          Do C = NOcc+1,NBF
            D = FindIndex(I,J,C)
            If(D <= NOcc) Cycle
            V_cdkl = ERI(C,D,K,L,dummy_flag=max(IRange,IRangeDirectLadders,IRangeQuadLadders))
            V_cdlk = ERI(C,D,L,K,dummy_flag=max(IRange,IRangeExchangeLadders,IRangeQuadLadders))
            If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
             Cycle
            End If
            J2_ijkl = J2_ijkl  + F12*(V_cdkl*T2abab(I,J,C) - V_cdlk*T2abba(I,J,C))
            J3_ijkl = J3_ijkl  + F12*(V_cdkl*T2abba(I,J,C) - V_cdlk*T2abab(I,J,C))
          End Do
! The intermediates are done, so contract with T2
          Do A = NOcc+1,NBF
            B = FindIndex(I,J,A)
            If(B <= NOcc) Cycle
            G2abab(I,J,A) = G2abab(I,J,A) + F12*(J2_ijkl*T2abab(K,L,A) + J3_ijkl*T2abba(K,L,A))
            G2abba(I,J,A) = G2abba(I,J,A) + F12*(J2_ijkl*T2abba(K,L,A) + J3_ijkl*T2abab(K,L,A))
          End Do
        End Do

! Add the particle-particle ladder
        Do A = NOcc+1,NBF
          B = FindIndex(I,J,A)
          If(B <= NOcc) Cycle
          Do C = NOcc+1,NBF
            D = FindIndex(I,J,C)
            If(D <= NOcc) Cycle
            V_cdab = ERI(C,D,A,B,dummy_flag=max(IRange,IRangeDirectLadders,IRangeLinLadders))
            V_cdba = ERI(C,D,B,A,dummy_flag=max(IRange,IRangeExchangeLadders,IRangeLinLadders))
            If(Min(V_cdab,V_cdba) < -80.0_pr) Then
             Print *, "Disallowed excitation not trapped!"
             Cycle
            End If
            G2abab(I,J,A) = G2abab(I,J,A) + F12*(V_cdab*T2abab(I,J,C) - V_cdba*T2abba(I,J,C))
            G2abba(I,J,A) = G2abba(I,J,A) + F12*(V_cdab*T2abba(I,J,C) - V_cdba*T2abab(I,J,C))
          End Do
        End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Ladders






      Subroutine GetG2Rings(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRange, &
                            IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF
      Real (Kind=pr), Intent(In)    :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abba(NOcc,NOcc,NOcc+1:NBF)
      Integer,        Intent(In)    :: IRange
      Integer,        Intent(In)    :: IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings
! Local variables
      Integer :: I, J, A, B, K, L, C, D
      Real (Kind=pr) :: J1_idal, J2_idal, J3_idal
      Real (Kind=pr) :: V_cjkb, V_cjbk, V_cika, V_ciak, V_cdkl, V_cdlk

!===========================================================!
!  Build the ring terms.  These are                         !
!   dG_ij^ab = t_{ik}^{ac} vbar^{kb}_{cj}                   !
!            + t_{kj}^{cb} vbar^{ak}_{ic}                   !
!            + t_{ik}^{ac} t_{lj}^{db} vbar^{kl}_{cd}       !
!===========================================================!

! These are the linear terms
      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle

! Do the t_ik^ac V_jc^bk contraction
        Do K = 1,NOcc
          C = FindIndex(I,K,A)
          If(C <= NOcc) Cycle
          V_cjkb = ERI(C,J,K,B,dummy_flag=max(IRange,IRangeDirectRings,IRangeLinRings))
          V_cjbk = ERI(C,J,B,K,dummy_flag=max(IRange,IRangeExchangeRings,IRangeLinRings))
          If(Min(V_cjkb, V_cjbk) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) + (V_cjkb - V_cjbk)*T2abab(I,K,A) +  V_cjkb*T2aaaa(I,K,A)
          G2abba(I,J,A) = G2abba(I,J,A) - V_cjbk*T2abba(I,K,A)
        End Do

! Do the t_jk^bc V_ic^ak contraction
        Do K = 1,NOcc
          C = FindIndex(J,K,B)
          If(C <= NOcc) Cycle
          V_cika = ERI(C,I,K,A,dummy_flag=max(IRange,IRangeDirectRings,IRangeLinRings))
          V_ciak = ERI(C,I,A,K,dummy_flag=max(IRange,IRangeExchangeRings,IRangeLinRings))
          If(Min(V_cika, V_ciak) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) + (V_cika - V_ciak)*T2abab(J,K,B) + V_cika*T2aaaa(J,K,B)
          G2abba(I,J,A) = G2abba(I,J,A) - V_ciak*T2abba(J,K,B)
        End Do
      End Do
      End Do
      End Do


! These are the quadratic terms.  They require intermediates.
! We'll have J1_idal, which contracts with T2aaaa(L,J,D,B) in G2aaaa(I,J,A,B) and with T2abab(L,J,D,B) in G2abab(I,J,A,B)
! We'll have J2_idal, which contracts with T2abab(L,J,D,B) in G2aaaa(I,J,A,B) and with T2aaaa(L,J,D,B) in G2abab(I,J,A,B
! We'll have J3_idal, which contracts with T2abba(L,J,D,B) in G2abba(I,J,A,B)
! This one is a true pain in the ass.
      Do I = 1,NOcc
      Do A = NOcc+1,NBF
        Do L = 1,NOcc
          J1_idal = Zero   ! Intermediate for T2aaaa(L,J,D,B) in G2aaaa(I,J,A,B) and T2abab(L,J,D,B) in G2abab(I,J,A,B)
          J2_idal = Zero   ! Intermediate for T2abab(L,J,D,B) in G2aaaa(I,J,A,B) and T2aaaa(L,J,D,B) in G2abab(I,J,A,B)
          J3_idal = Zero   ! Intermediate for T2abba(L,J,D,B) in G2abba(I,J,A,B)
          Do K = 1,NOcc
            C = FindIndex(I,K,A)
            If(C <= NOcc) Cycle
            D = FindIndex(K,L,C)
            If(D <= NOcc) Cycle
            V_cdkl = ERI(C,D,K,L,dummy_flag=max(IRange,IRangeDirectRings,IRangeQuadRings))
            V_cdlk = ERI(C,D,L,K,dummy_flag=max(IRange,IRangeExchangeRings,IRangeQuadRings))
            If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
              Cycle
            End If
            J1_idal = J1_idal + T2aaaa(I,K,A)*(V_cdkl - V_cdlk) + T2abab(I,K,A)*V_cdkl
            J2_idal = J2_idal + T2abab(I,K,A)*(V_cdkl - V_cdlk) + T2aaaa(I,K,A)*V_cdkl
            J3_idal = J3_idal - T2abba(I,K,A)*V_cdlk
          End Do

! The intermediates are done, so contract with T2
! In principle, we've got the right D here.
          Do J = 1,NOcc
            B = FindIndex(I,J,A)
            If(B <= NOcc) Cycle
            D = FindIndex(J,L,B)
            If(D <= NOcc) Cycle
            G2abab(I,J,A) = G2abab(I,J,A) + J2_idal*T2aaaa(L,J,D) + J1_idal*T2abab(L,J,D)
            G2abba(I,J,A) = G2abba(I,J,A) + J3_idal*T2abba(L,J,D)
          End Do
        End Do
      End Do
      End Do

      Return
      End Subroutine GetG2Rings






      Subroutine GetG2XRings(G2abab,G2abba,T2aaaa,T2abab,T2abba,NOcc,NBF,IRange, &
                             IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings)
      Implicit None
      Integer,        Intent(In)    :: NOcc, NBF
      Real (Kind=pr), Intent(In)    :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(In)    :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(InOut) :: G2abba(NOcc,NOcc,NOcc+1:NBF)
      Integer,        Intent(In)    :: IRange
      Integer,        Intent(In)    :: IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings
! Local variables
      Integer :: I, J, A, B, K, L, C, D
      Real (Kind=pr) :: J1_idlb, J2_idlb, J3_idlb
      Real (Kind=pr) :: V_jcka, V_ickb, V_jcak, V_icbk, V_cdkl, V_cdlk

!===========================================================!
!  Build the crossed-ring terms.  These are                 !
!   dG_ij^ab = - t_{ik}^{cb} vbar^{ka}_{jc}                 !
!              - t_{jk}^{ca} vbar^{kb}_{ic}                 !
!              - t_{ik}^{cb} t_{jl}^{da} vbar^{kl}_{cd}     !
!  Most of the signs drop out on using vbar!                !
!===========================================================!

! These are the linear terms
      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle

! Do the t_ik^cb V_jc^ka contraction
        Do K = 1,NOcc
          C = FindIndex(K,I,B)
          If(C <= NOcc) Cycle
          V_jcka = ERI(J,C,K,A,max(IRange,IRangeExchangeRings,IRangeLinRings))
          V_jcak = ERI(J,C,A,K,max(IRange,IRangeDirectRings,IRangeLinRings))
          If(Min(V_jcka,V_jcak) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) - V_jcka*T2abab(K,I,B)
          G2abba(I,J,A) = G2abba(I,J,A) + (V_jcak - V_jcka)*T2abba(K,I,B)  +  V_jcak*T2aaaa(K,I,B)
        End Do

! Do the t_jk^ca v_ic^kb contraction
        Do K = 1,NOcc
          C = FindIndex(K,J,A)       ! should return the same as GetIndex(I,K,B), but check!
          If(C <= NOcc) Cycle
          V_ickb = ERI(I,C,K,B,max(IRange,IRangeExchangeRings,IRangeLinRings))
          V_icbk = ERI(I,C,B,K,max(IRange,IRangeDirectRings,IRangeLinRings))
          If(Min(V_ickb,V_icbk) < -80.0_pr) Then
            Print *, "Disallowed excitation not trapped!"
            Cycle
          End If
          G2abab(I,J,A) = G2abab(I,J,A) - V_ickb*T2abab(K,J,A)
          G2abba(I,J,A) = G2abba(I,J,A) + (V_icbk - V_ickb)*T2abba(K,J,A) + V_icbk*T2aaaa(K,J,A)
        End Do
      End Do
      End Do
      End Do


! These are the quadratic terms.  They require intermediates.
! We'll have J1_idlb, which contracts with T2aaaa(L,J,A,D) in G2aaaa(I,J,A,B) and with T2abba(L,J,A,D) in G2abba(I,J,A,B)
! We'll have J2_idlb, which contracts with T2abba(L,J,A,D) in G2aaaa(I,J,A,B) and with T2aaaa(L,J,A,D) in G2abba(I,J,A,B
! We'll have J3_idlb, which contracts with T2abab(L,J,A,D) in G2abab(I,J,A,B)
      Do I = 1,NOcc
      Do B = NOcc+1,NBF
        Do L = 1,NOcc
          J1_idlb = Zero   ! Intermediate for T2aaaa(L,J,A,D) in G2aaaa(I,J,A,B) and T2abba(L,J,A,D) in G2abba(I,J,A,B)
          J2_idlb = Zero   ! Intermediate for T2abba(L,J,A,D) in G2aaaa(I,J,A,B) and T2aaaa(L,J,A,D) in G2abba(I,J,A,B
          J3_idlb = Zero   ! Intermediate for T2abab(L,J,D,B) in G2abab(I,J,A,B)
          Do K = 1,NOcc
            C = FindIndex(I,K,B)
            If(C <= NOcc) Cycle
            D = FindIndex(K,L,C)
            If(D <= NOcc) Cycle
            V_cdkl = ERI(C,D,K,L,max(IRange,IRangeDirectRings,IRangeQuadRings))
            V_cdlk = ERI(C,D,L,K,max(IRange,IRangeExchangeRings,IRangeQuadRings))
            If(Min(V_cdkl,V_cdlk) < -80.0_pr) Then
              Print *, "Disallowed excitation not trapped!"
              Cycle
            End If
            J1_idlb = J1_idlb + T2aaaa(K,I,B)*(V_cdkl - V_cdlk) + T2abba(K,I,B)*V_cdkl
            J2_idlb = J2_idlb + T2abba(K,I,B)*(V_cdkl - V_cdlk) + T2aaaa(K,I,B)*V_cdkl
            J3_idlb = J3_idlb + T2abab(K,I,B)*V_cdlk
          End Do
! The intermediates are done, so contract with T2
          Do J = 1,NOcc
            A = FindIndex(J,I,B)
            If(A <= NOcc) Cycle
            D = FindIndex(L,J,A)
            If(D <= NOcc) Cycle
            G2abba(I,J,A) = G2abba(I,J,A) - J1_idlb*T2abba(L,J,A) - J2_idlb*T2aaaa(L,J,A)
            G2abab(I,J,A) = G2abab(I,J,A) + J3_idlb*T2abab(L,J,A)
          End Do
        End Do
      End Do
      End Do

      Return
      End Subroutine GetG2XRings

   End Module CCD

