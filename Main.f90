
      Program TermedCC
      Use Precision
      Use IO
      Use MP2
      Use CCD
      Implicit None
! Dimensioning variables
      Integer :: NOcc, NAO, NumRSPoints, iRSPoint
      Logical :: DoRing, DoXRing, DoLadder, DoMosaic
      Integer :: IRangeRing, IRangeXRing, IRangeLadder, IRangeMosaic
      Integer :: IRangeDriverDirect, IRangeDriverExchange, IRangeEnergy
      Integer :: IRangeLinRings, IRangeQuadRings, IRangeDirectRings, IRangeExchangeRings  
      Integer :: IRangeLinLadders, IRangeQuadLadders, IRangeDirectLadders, IRangeExchangeLadders
      Real (Kind=pr) :: rSMin, rSMax, rS
! Correlated stuff
      Real (Kind=pr) :: ECorr
      Real (Kind=pr), Allocatable :: T2aaaa(:,:,:), T2abab(:,:,:), T2abba(:,:,:)
      Real (Kind=pr), Allocatable :: X2aaaa(:,:,:), X2abab(:,:,:), X2abba(:,:,:)
! Error checking variables
      Integer, Parameter :: NAlloc = 6
      Integer :: IAlloc(NAlloc)
      Logical, Parameter :: T = .true., F=.false.

!==========================================!
!  This code implements RHF-based CCD.     !
!  We give it the option of keeping the    !
!  ring, ladder, crossed-ring, and mosaic  !
!  terms each on a case-by-case basis.     !
!------------------------------------------!
!  The first thing we must do is read the  !
!  basic information for the calculation.  !
!==========================================!

      Call ReadInput(NOcc,NAO,rSMin,rSMax,NumRSPoints,       &
                     DoRing,DoXRing,DoLadder,DoMosaic,       &
                     IRangeRing,IRangeXRing,IRangeLadder,IRangeMosaic, &
                     IRangeDriverDirect,IRangeDriverExchange,IRangeEnergy,                       &
                     IRangeLinRings,IRangeQuadRings,IRangeDirectRings,IRangeExchangeRings, &
                     IRangeLinLadders,IRangeQuadLadders,IRangeDirectLadders,IRangeExchangeLadders)


!==========================================!
!  Now we can allocate the memory and go!  !
!==========================================!

      IAlloc = 0
      Allocate(T2aaaa(NOcc,NOcc,NOcc+1:NAO),  Stat=IAlloc(1))
      Allocate(T2abab(NOcc,NOcc,NOcc+1:NAO),  Stat=IAlloc(2))
      Allocate(T2abba(NOcc,NOcc,NOcc+1:NAO),  Stat=IAlloc(3))
      Allocate(X2aaaa(NOcc,NOcc,NOcc+1:NAO),  Stat=IAlloc(4))
      Allocate(X2abab(NOcc,NOcc,NOcc+1:NAO),  Stat=IAlloc(5))
      Allocate(X2abba(NOcc,NOcc,NOcc+1:NAO),  Stat=IAlloc(6))
      If(Any(IAlloc /= 0)) Stop "Could not allocate in main"
      Open(7,File='Output',Status="Replace")
      Close(7)


!===============================!
!  Loop over rS values and go.  !
!===============================!

      Do iRSPoint = 1,NumRSPoints
        If (NumRSPoints.eq.1) Then
          rS=rSmin
        Else
          rS = rSMin + (iRSPoint-1)*(rSMax-rSMin)/(NumRSPoints-1)
        End If
        Call change_rs(rS)
        Call DrvMBPT(Eigen,X2aaaa,X2abab,X2abba,NOcc,NAO,EHF,ECorr)
        If(iRSPoint == 1) Then
          T2aaaa = X2aaaa
          T2abab = X2abab
          T2abba = X2abba
        End If
        Call DrvCCD(Eigen,T2aaaa,T2abab,T2abba,NOcc,NAO,EHF,ECorr,  &
                    DoRing,DoXRing,DoLadder,DoMosaic, &
                    IRangeRing,IRangeXRing,IRangeLadder,IRangeMosaic, &
                    IRangeDriverDirect,IRangeDriverExchange,IRangeEnergy,                       &
                    IRangeLinRings,IRangeQuadRings,IRangeDirectRings,IRangeExchangeRings, &
                    IRangeLinLadders,IRangeQuadLadders,IRangeDirectLadders,IRangeExchangeLadders)
      End Do


!==============================================!
!  Lastly, deallocate memory and exit safely.  !
!==============================================!

      DeAllocate(T2aaaa,    Stat=IAlloc(1))
      DeAllocate(T2abab,    Stat=IAlloc(2))
      DeAllocate(T2abba,    Stat=IAlloc(3))
      DeAllocate(X2aaaa,    Stat=IAlloc(4))
      DeAllocate(X2abab,    Stat=IAlloc(5))
      DeAllocate(X2abba,    Stat=IAlloc(6))
      If(Any(IAlloc /= 0)) Stop "Could not deallocate in main"

      Stop
      End Program TermedCC

