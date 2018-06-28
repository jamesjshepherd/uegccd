module heg 

use Precision
use Constants

implicit none

integer :: ecut!=9       ! In scaled units n^2, maximum is cell^2
integer :: Nel!=54       ! Number of electrons
real(pr) :: rs!=1.0_pr  ! rs
   Real (Kind=pr) :: EHF

type basis_set
    integer :: sorted_index
    real(pr) :: k2
    integer :: n2
    integer :: n(3)
    integer :: ms
end type basis_set

real(pr) :: L, om, madelung
real(pr) :: fermi_wavevector, density, screening_distance
integer :: cell
integer :: nBasis
real(pr), allocatable :: eigen(:) 
real(pr), allocatable :: kvec(:,:)
integer, allocatable :: nvec(:,:)
integer, allocatable :: kPointToBasisFn(:,:,:)

private
public :: eigen, ERI, init_HEG, FindIndex, nBasis, EHF, init_HEG_dummy, change_rs, findtol

contains

    subroutine init_HEG(ecut_tmp,Nel_tmp,rs_tmp)
        
        integer, intent(in) :: nEl_tmp, ecut_tmp
        real(pr), intent(in) :: rs_tmp

        call set_up_data(ecut_tmp,Nel_tmp,rs_tmp)
        call print_shells(10) 
        call setup
        call setup_basis

    end subroutine

    
    
    
    
    
    subroutine set_up_data(ecut_tmp,Nel_tmp,rs_tmp)

        integer :: nEl_tmp, ecut_tmp
        real(pr) :: rs_tmp

        Nel=Nel_tmp
        rs=rs_tmp
        ecut=ecut_tmp

    end subroutine






    subroutine print_shells(kc) 
    ! kc is an interger input
    ! genuinely kc in integer coordinates
    ! corresponds to "cell" keyword

        integer :: l1,l2,l3
        integer :: kc
        integer :: bins
        integer :: energy
        integer, allocatable :: histogram(:)
        integer :: culm

        bins=kc**2
        allocate(histogram(0:bins))
        histogram = 0

        do l1=-kc,kc
            do l2=-kc,kc
                do l3=-kc,kc
                    energy=l1**2+l2**2+l3**2
                    if (energy.lt.kc**2) then
                        histogram(energy)=histogram(energy)+1
                    endif
                enddo
            enddo
        enddo
        culm=0
        write(60,*) ""
        write(60,*) "--------------------"
        write(60,*) "Basis sets available"
        write(60,*) "--------------------"
        write(60,*) "# ecut k-points M"

        do l1=0,kc**2
            if (histogram(l1).ne.0) then
                culm=culm+histogram(l1)
                write(60,*) l1, culm, culm*2
            endif
        enddo

    end subroutine






    subroutine setup

        L=rs*(4.0_pr/3.0_pr*pi*Nel)**(1.0_pr/3.0_pr)
        om=L**3
        madelung=calc_madelung()
        ! The condition on cell is: cell > sqrt(ecut)
        cell=int(sqrt(real(ecut,dp)))+1 
        write(60,*) "Madelung constant calculated as", madelung
        density=Nel/L**3.0_pr
        fermi_wavevector=(3.0_pr*pi**2.0_pr*density)**(1.0_pr/3.0_pr)
        screening_distance=(4.0_dp/pi*fermi_wavevector)**0.5_dp

    end subroutine 






    function calc_madelung()
        
        real(pr) :: calc_madelung
        real(pr) :: kappa
        integer :: l1,l2,l3,l4
        integer :: n2
        real(pr) :: k2,ek2,recipsum2
        real(pr) :: modr,er2,realsum2
        real(pr) :: term2,term4

        write(60,*) ""
        write(60,*) "-----------------------------------------------------------"
        write(60,*) "Calculating Madelung Constant - Fraser et al. PRB 53 4 1814"
        write(60,*) "-----------------------------------------------------------"
        kappa=1.0_pr/om**(1.0_pr/3.0_pr) ! there are more optimal values, but this is fine for here

        term2=-pi/(kappa**2.0_pr*L**3.0_pr)
        write(60,*) term2, "term2"
        term4=-2.0_pr*kappa/sqrt(pi)
        write(60,*) term4, "term4"

        recipsum2=0.0_pr
        do l4=1,10
            recipsum2=0.0_pr
            do l1=-l4,l4
                do l2=-l4,l4
                    do l3=-l4,l4
                        n2=l1**2+l2**2+l3**2
                        k2=(1.0_pr/L**2.0_pr)*(l1**2+l2**2+l3**2)
                        ek2=(1.0_pr/L**3.0_pr)*(1.0_pr/(pi*k2))*exp(-pi**2.0_pr*k2/kappa**2.0_pr)
                        if (n2.ne.0) then
                            !write(60,*) k2,ek2 ! for testing
                            recipsum2=recipsum2+ek2
                        endif
                    enddo
                enddo
            enddo
            write(60,*) l4,recipsum2
        enddo
        write(60,*) "reciprocal space", recipsum2
        
        realsum2=0.0_pr
        do l4=1,10
            realsum2=0.0_pr
            do l1=-l4,l4
                do l2=-l4,l4
                    do l3=-l4,l4
                        modr=L*sqrt(real((l1**2+l2**2+l3**2),dp))
                        if (modr.ne.0.0_pr) then
                            er2=erfc(kappa*modr)/modr
                            realsum2=realsum2+er2
                            !write(60,*) modr*0.5_pr,er2 ! for testing
                        endif
                    enddo
                enddo
            enddo
            write(60,*) l4,realsum2
        enddo
        write(60,*) "real space", realsum2
        calc_madelung=realsum2+recipsum2+term2+term4
        !calc_madelung=0.0_pr!realsum2+recipsum2+term2+term4
        madelung=calc_madelung

    end function

    subroutine setup_basis

        type (basis_set), allocatable :: G1(:)
        type (basis_set) :: G1_temp
        logical :: swapped
        integer :: l1, l2, l3
        integer :: M
        real(pr) :: qvec(3)

        write(60,*) ""
        write(60,*) "-----------------------------------------------------------"
        write(60,*) "Setting up HEG basis set"
        write(60,*) "-----------------------------------------------------------"
        allocate(G1((2*cell+1)**3))
        write(60,*) "Memory allocated for", (2*cell+1)**3, "basis functions"

        ! Generates basis set: all plane waves which have an less than the energy 
        ! cutoff ecut where E=(1/2)(n^2) i.e. scaled units

        ! Plane wave bais set:
        ! \psi_j(r,\sigma) = \sqrt(1/\Omega)*exp(i*k_i.r) \delta_{\sigma_i,\sigma}
        ! where real-space cell volume \Omega = L**3 and 
        ! the reciprocal lattice vectors are k = ((2*pi)/L)(n,m,l) 
        ! n,m,l=0,+/-1,+/-2,... are integers

        M=0
        do l1=-cell,cell
            do l2=-cell,cell
                do l3=-cell,cell
                    M=M+1
                    G1(M)%n(1)=l1
                    G1(M)%n(2)=l2
                    G1(M)%n(3)=l3
                    G1(M)%n2=l1**2+l2**2+l3**2
                    G1(M)%k2=(l1**2+l2**2+l3**2)*(2*pi/L)**2
                    if (G1(M)%n2 .gt. ecut) then 
                        M=M-1 ! By construction G1(M+1) exists but is never used (caution!)
                    endif
                enddo
            enddo
        enddo

        nBasis=M

        write(60,*) "Actual basis set size", nBasis

        ! Sorting basis set according to ascending energy
        do
            swapped=.false.
            do l1=2,nBasis
                if (G1(l1)%n2.lt.G1(l1-1)%n2) then
                    G1_temp=G1(l1-1)
                    G1(l1-1)=G1(l1)
                    G1(l1)=G1_temp
                    swapped=.true.
                endif
            enddo
            if (.not.swapped) exit
        enddo


        ! Initiatlise arrays
        allocate(kvec(M,3))
        allocate(eigen(M))
        allocate(nvec(M,3))
        kvec=0.0_pr
        eigen=0.0_pr
        nvec=0

        ! Form the array of k vectors (scaled)
        ! and n vectors (integers)
        do l1=1,nbasis 
            kvec(l1,1)=G1(l1)%n(1)*(2.0_pr*pi/L)
            kvec(l1,2)=G1(l1)%n(2)*(2.0_pr*pi/L)
            kvec(l1,3)=G1(l1)%n(3)*(2.0_pr*pi/L)
            nvec(l1,1)=G1(l1)%n(1)
            nvec(l1,2)=G1(l1)%n(2)
            nvec(l1,3)=G1(l1)%n(3)
        enddo

        ! Form eigenvector array
        do l1=1,nbasis
            eigen(l1)=0.0_pr
            eigen(l1)=0.5_pr*dot_product(kvec(l1,:),kvec(l1,:)) ! kinetic
            do l2=1,nEl/2
                if (l2.ne.l1) then
                    ! find "q" momentum transfer
                    qvec=kvec(l1,:)-kvec(l2,:)
                    eigen(l1)=eigen(l1)-(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))
                else 
                    eigen(l1)=eigen(l1)+madelung ! self-interaction, madelung is a negative constant
                endif
            enddo
            write(60,*) l1, nvec(l1,1), nvec(l1,2), nvec(l1,3), eigen(l1)
        enddo
      
        ! Print out E_HF
        call print_HF

        ! Create a look-up table from the k-points to the orbital indices
        allocate(kPointToBasisFn(-cell:cell,-cell:cell,-cell:cell))
        ! write(60,*) "Memory allocated for", (2*cell+1)**3, "k point look-up table"
        kPointToBasisFn=0 ! there are some elements which aren't in the basis, these will
                          ! be returned as zero
        do l1=1,nBasis
            kPointToBasisFn(  G1(l1)%n(1), &
                              G1(l1)%n(2), &
                              G1(l1)%n(3)     ) = l1
        enddo

        ! G1 is now an unnecessary array
        deallocate(G1)

    end subroutine setup_basis

    ! This function takes integer inputs as orbital indices and finds the fourth 
    ! momentum-allowed index and returns the value in the function
    !
    ! Returns FindIndex=0 if there isn't a momentum allowed orbital in the WHOLE BASIS
    !
    ! The momentum conservation law used here is that for an excitation 
    ! ij -> ab
    ! k_i+k_j = k_a+k_b
    ! i.e. the total momentum is the same before and after the excitation
    !
    ! Inputs: m,n,l - indecies of orbitals
    ! Global data requirements: kPointToBasisFn, nvec, ecut, nBasis
    !
    function FindIndex(m,n,l)

        integer :: FindIndex
        integer :: m,n,l
        integer :: kTarget(3)
        integer :: eTarget
        logical :: InBasisSet

        kTarget=nvec(m,:)+nvec(n,:)-nvec(l,:)

        ! Is this in the basis set at all?
       InBasisSet=.true.
        ! First, check range
        if (abs(kTarget(1)).gt.cell) InBasisSet=.false. 
        if (abs(kTarget(2)).gt.cell) InBasisSet=.false. 
        if (abs(kTarget(3)).gt.cell) InBasisSet=.false. 
        ! Second, check energy
        eTarget=dot_product(kTarget,kTarget)
        if (eTarget.gt.ecut) InBasisSet=.false.
        ! If it is out then return zero
        if (.not.InBasisSet) then
            FindIndex=0
            return
        endif

        ! Now it's safe to use the look-up table
        FindIndex=kPointToBasisFn(kTarget(1),kTarget(2),kTarget(3))

        if ((FindIndex.eq.0) .or. (FindIndex.gt.nBasis))  then
            stop "Error in FindIndex" ! This is now an error since it should 
            ! have been filtered out before it hit the lookup table. We check 
            ! this because the table might have array bound errors
        endif

    end function

    ! As FindIndex(m,n,l) but only returns non-zero index if it's also in the 
    ! occupied manifold
    !
    function FindIndexOcc(m,n,l)

        integer :: FindIndexOcc
        integer :: m,n,l

        FindIndexOcc=FindIndex(m,n,l)
        if (FindIndexOcc.gt.nEl/2) then 
            FindIndexOcc=0
        endif

    end function

    ! As FindIndex(m,n,l) but only returns non-zero index if it's also in the 
    ! virtual manifold
    !
    function FindIndexVirt(m,n,l)

        integer :: FindIndexVirt
        integer :: m,n,l

        FindIndexVirt=FindIndex(m,n,l)
        if (.not.(FindIndexVirt.le.nEl/2)) then ! i.e. if it's not in the OCCUPIED manifold
            FindIndexVirt=0
        endif

    end function

    ! This routine, principally a debug routine, checks for momentum symmetry of four indices
    ! 
    function check_momentum_symmetry(i,j,a,b)

        logical :: check_momentum_symmetry
        integer :: i, j, a, b
        integer :: g1(3), g2(3)

        check_momentum_symmetry=.false.
        g1=-nvec(i,:)+nvec(a,:) ! i.e. g1 goes from i to a
        g2=-nvec(j,:)+nvec(b,:) ! i.e. g2 goes from j to b
        if (g1(1).eq.-g2(1).and.g1(2).eq.-g2(2).and.g1(3).eq.-g2(3)) check_momentum_symmetry=.true.

    end function

    ! This takes an INTEGER VECTOR input and returns the electron repulsion integral
    !
    function ERI1(qvec)

        real(pr) :: ERI1
        integer :: qvec(3)

        ERI1=(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))

    end function

    ! This takes FOUR INTEGER ORBITAL INDICES and returns the electron repulsion integral
    ! for an excitation ij->ab
    ! Strangely, ERI is returned simply as a negative number if it's momentum disallowed
    ! (when the integral is actually zero)
    !
!    function ERI(a,b,i,j)
!
!        integer :: i,j,a,b
!        real(pr) :: ERI
!        integer :: qvec(3)
!
!        if (min(a,b,i,j) <= 0) then
!            stop "ERROR: this is not an allowed excitation"
!        endif
!        if (.not.check_momentum_symmetry(i,j,a,b)) then
!            ERI=-100.0_pr
!!           write(60,*) "WARNING: ERI called as zero, returned as -100!"
!            return
!        endif
!        qvec=-nvec(i,:)+nvec(a,:)
!        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
!            ERI=-madelung
!            !stop "ERROR: zero momentum error"
!        else
!            !ERI=(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))
!            ERI= One/(L*Pi*Dot_Product(qVec,qVec))
!        endif
!
!    end function

    
    subroutine init_HEG_dummy(ecut_tmp,Nel_tmp)
    ! This is a dummy routine ahead of an rs scan version of this code
    ! so rs is initalised as one (because scale factors are simpler)
        
        integer, intent(in) :: nEl_tmp, ecut_tmp

        call init_HEG(ecut_tmp,Nel_tmp,1.0_pr)

    end subroutine

    subroutine change_rs(new_rs)

        real(pr), intent(in) :: new_rs
        real(pr) :: old_rs, rescale
        real(pr) :: qvec(3)
        integer :: l1, l2

        old_rs=rs
        rs=new_rs
        rescale=new_rs/old_rs 
        ! imagine old_rs as 1 and everything makes sense...
        L=L*rescale
        om=om*rescale**3.0_pr
        kvec=kvec*(1.0_pr/rescale)
        madelung=calc_madelung()
        density=Nel/L**3.0_pr
        fermi_wavevector=(3.0_pr*pi**2.0_pr*density)**(1.0_pr/3.0_pr)
        screening_distance=(4.0_dp/pi*fermi_wavevector)**0.5_dp


        ! Form eigenvector array
        eigen=0.0_pr
        do l1=1,nbasis
            eigen(l1)=0.0_pr
            eigen(l1)=0.5_pr*dot_product(kvec(l1,:),kvec(l1,:)) ! kinetic
            do l2=1,nEl/2
                if (l2.ne.l1) then
                    ! find "q" momentum transfer
                    qvec=kvec(l1,:)-kvec(l2,:)
                    eigen(l1)=eigen(l1)-(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))
                else 
                    eigen(l1)=eigen(l1)+madelung ! self-interaction, madelung is a negative constant
                endif
            enddo
            write(60,*) l1, nvec(l1,1), nvec(l1,2), nvec(l1,3), eigen(l1)
        enddo
        
        call print_HF

    end subroutine

    subroutine print_HF
        
        real(pr) :: e_hf
        real(pr) :: e_ex
        integer :: l1, l2
        real(pr) :: qvec(3)

        ! Print out E_HF
        e_hf=0.0_pr
        e_ex=0.0_pr
        do l1=1,nEl/2
            e_hf=e_hf+0.5_pr*dot_product(kvec(l1,:),kvec(l1,:))
            do l2=1,nEl/2
                if (l2.ne.l1) then
                    qvec=kvec(l1,:)-kvec(l2,:)
                    e_hf=e_hf-0.5_pr*(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))
                    e_ex=e_ex-0.5_pr*(1/om)*(4.0_pr*pi/dot_product(qvec,qvec))
                endif
            enddo
        enddo
        e_hf=e_hf*2.0_pr
        e_ex=e_ex*2.0_pr
        write(60,*) "L:", L
        write(60,*) "om:", om
        write(60,*) "cell:", cell
        write(60,*) "HF energy:", e_hf/nEl+madelung/2.0
        write(60,*) "Exchange energy:", e_ex/nEl+madelung/2.0
        write(60,*) "Madelung:", madelung
        EHF = e_hf + madelung*F12*nEL

    end subroutine

    ! This takes FOUR INTEGER ORBITAL INDICES and returns the electron repulsion integral
    ! for an excitation ij->ab
    ! Strangely, ERI is returned simply as a negative number if it's momentum disallowed
    ! (when the integral is actually zero)
    ! NEW: now has a flag argument for calculation of other kernels
    !
    function ERI(a,b,i,j,dummy_flag)

        integer :: i,j,a,b
        integer, intent(in), optional :: dummy_flag ! This is 1=short-range or 2=long-range
        integer :: flag ! This is 0=normal, 1=short-range or 2=long-range
        real(pr) :: ERI
        integer :: qvec(3)
        real(pr) :: qvec2, qvec2_scaled

        if (.not.present(dummy_flag)) then
            flag=0
        else
            flag=dummy_flag
        endif

        if (min(a,b,i,j) <= 0) then
            stop "ERROR: this is not an allowed excitation"
        endif
        if (.not.check_momentum_symmetry(i,j,a,b)) then
            ERI=-100.0_pr
            return
        endif
        qvec=-nvec(i,:)+nvec(a,:)
        if ((qvec(1).eq.0).and.(qvec(2).eq.0).and.(qvec(3).eq.0)) then
            select case (flag)
                case (0) ! normal
                    ERI=-madelung
                case (1) ! short-range
                    ERI=(1.0_dp/L**3.0_dp)*4.0_pr*pi/screening_distance**2.0_pr
                case (2) ! long-range = normal-short range
                    ERI=-madelung-(1.0_dp/L**3.0_dp)*4.0_pr*pi/screening_distance**2.0_pr
                case (3) ! completely screened, zero
                    ERI=0.0_pr
                case default
                    stop "ERI called incorrectly"
            end select
        else
            qvec2=dot_product(qvec,qvec)
            qvec2_scaled=qvec2*(2.0_pr*pi/L)**2.0_pr
            select case (flag)
                case (0) ! normal
                    ERI= One/(L*Pi*qvec2)
                case (1) ! short-range
                    ERI= (1.0_dp/L**3.0_dp)*4.0_dp*pi/(screening_distance**2.0_dp+qvec2_scaled)
                case (2) ! long-range = normal-short range
                    ERI= (1.0_dp/L**3.0_dp)*4.0_dp*pi/(screening_distance**2.0_dp+qvec2_scaled) &
                        *screening_distance**2.0_dp/qvec2_scaled 
                case (3) ! completely screened, zero
                    ERI=0.0_pr
                case default
            end select
        endif
        if (ERI.lt.0.0_pr) then
            stop "ERI less than zero"
        endif

    end function

    subroutine findtol(tolmax,failratio)

        real(pr) :: tolmax,failratio

        TolMax = 1.0E-8_pr
        FailRatio = 100
        if (rs.gt.10.1_pr) TolMax = 1.0E-5_pr
        if (rs.lt.0.95_pr) FailRatio = 1E6

    end subroutine
    
end module 
