
   Module MP2
   Use Precision
   Use Constants
   Use HEG

   Contains

      Subroutine DrvMBPT(Fock,T2aaaa,T2abab,T2abba,NOcc,NBF,ESCF,ECorr)
      Implicit None
      Integer,        Intent(In)  :: NOcc, NBF
      Real (Kind=pr), Intent(In)  :: ESCF
      Real (Kind=pr), Intent(In)  :: Fock(NBF)
      Real (Kind=pr), Intent(Out) :: T2aaaa(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(Out) :: T2abab(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(Out) :: T2abba(NOcc,NOcc,NOcc+1:NBF)
      Real (Kind=pr), Intent(Out) :: ECorr
      Real (Kind=pr) :: V_ijab, V_ijba, Denom
      Integer :: I, J, A, B

!==========================================!
!  This does the MP2 calculation.  Glory!  !
!==========================================!

      Write(6,*) "Doing MBPT2..."
      T2abab = Zero
      T2abba = Zero
      T2aaaa = Zero
      ECorr = Zero
      Do I = 1,NOcc
      Do J = 1,NOcc
      Do A = NOcc+1,NBF
        B = FindIndex(I,J,A)
        If(B <= NOcc) Cycle
        V_ijab = ERI(I,J,A,B)
        V_ijba = ERI(I,J,B,A)
        If(Min(V_ijab,V_ijba) < -80.0_pr) Then
          Print *, "Disallowed excitation not trapped!"
          Cycle
        End If
        Denom = Fock(I) + Fock(J) - Fock(A) - Fock(B)
        T2abab(I,J,A) =  V_ijab/Denom
        T2abba(I,J,A) = -V_ijba/Denom
        T2aaaa(I,J,A) =  T2abab(I,J,A) + T2abba(I,J,A)
        ECorr = ECorr                                                &
              + F12*T2aaaa(I,J,A)*(ERI(I,J,A,B) - ERI(I,J,B,A))      &
              + F12*T2abab(I,J,A)*ERI(I,J,A,B)                       &
              - F12*T2abba(I,J,A)*ERI(I,J,B,A)
      End Do
      End Do
      End Do


! Write to Output and we're done!
      Open(7,File='Output',Position="Append")
      Write(7,1000)
      Write(7,1005) ESCF
      Write(7,1040)
      Write(7,1030)
      Write(7,1040)
      Do I = 1,NOcc
        Write(7,1050) Fock(I)
      End Do
      Write(7,1060)
      Do A = NOcc+1,NBF
        Write(7,1050) Fock(A)
      End Do
      Write(7,1000)
      Write(7,1010) ECorr
      Write(7,1020) ECorr + ESCF
      Write(7,1000)
      Close(7)


1000  Format(14x,'**************************************************')
1005  Format(14x,'* E(SCF) is ',12x,F18.12,' a.u.',2x,'*')
1010  Format(14x,'* E(2) is ',14x,F18.12,' a.u.',2x,'*')
1020  Format(14x,'* E(SCF) + E(2) is ',5x,F18.12,' a.u.',2x,'*')
1030  Format(14X,'*                RHF Eigenvalues:                *')
1040  Format(14x,'*------------------------------------------------*')
1050  Format(14x,'*',12x,ES24.16,12x,'*')
1060  Format(14x,'*',12x,24("-"),12x,'*')

      Return
      End Subroutine DrvMBPT

   End Module MP2

