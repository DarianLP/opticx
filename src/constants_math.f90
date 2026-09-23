module constants_math
  implicit none
  real(8), parameter :: pi=3.14159265358979323846d0
  real(8), parameter :: dk=1.0d-6
  
  contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine percentage_index(kacum,ktotal,kmoment)
    implicit none
    integer :: kacum,ktotal,kmoment
    integer :: npercentage,nrest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    npercentage=int(dble(kacum)/dble(ktotal)*100.0d0)
    nrest=mod(npercentage,10)
    !write(*,*) 'Percentage of k-points read:',npercentage,' %'
    if (nrest.eq.0) then
      if (kmoment.ne.npercentage) then
        write(*,*) '   Percentage of loop:',npercentage,' %'
      end if
      kmoment=npercentage
    end if   
  end subroutine percentage_index
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine crossproduct(ax,ay,az,bx,by,bz,cx,cy,cz)
    implicit none
    real(8) :: ax,ay,az,bx,by,bz,cx,cy,cz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cx=ay*bz-az*by
    cy=az*bx-ax*bz
    cz=ax*by-ay*bx     
  end subroutine crossproduct


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         diagoz
!   INPUTS:       h matrix to diagonalize
!                 n dimension of h
!   OUTPUTS:      w; eigenvalues of h
!                 h;  gives eigenvectors by columns as output
!   DESCRIPTION:  this subroutine uses Lapack libraries to diagonalize
!                 an hermitian complex matrix.
!   
!     Juan Jose Esteve-Paredes                28.11.2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine diagoz(n,w,h)
    !finding the eigenvalues of a complex matrix using LAPACK
    implicit real*8 (a-h,o-z)
    !declarations, notice double precision
    integer n,INFO,LWORK
    dimension w(n)
    dimension RWORK(3*n-2)
    dimension h(n,n)

    real(8) w
    complex*16 h
    complex*16 WORK(2*n)
    character*1 JOBZ,UPLO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !find the solution using the LAPACK routine ZGEEV
    JOBZ='V'
    UPLO='U'
    LWORK=2*n
            
    call zheev(JOBZ, UPLO, n, h, n, w, WORK, LWORK, RWORK, INFO)
  end

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         diagoz_gen
!   INPUTS:       h: matrix to diagonalize
!                 s: overlap matrix 
!                 n: dimension of h and s
!   OUTPUTS:      w; eigenvalues of h
!                 h;  gives eigenvectors by columns as output
!   DESCRIPTION:  this subroutine solves the generalized eigenvalue problem H*v = e*S*v using 
!                 LAPACK zhegv
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine diagoz_gen(n,w,h,s)
    implicit none 
    integer n,INFO,LWORK,ITYPE
    dimension w(n)
    dimension RWORK(3*n-2)
    dimension h(n,n)
    dimension s(n,n)

    real(8) w
    real(8) RWORK
    complex*16 h,s
    complex*16 WORK(2*n)
    character*1 JOBZ,UPLO

    ! diagnostic variables
    integer :: info_tmp, i, j
    real(8) :: s_eigs(n), rwork_tmp(3*n-2)
    complex*16 :: s_copy(n,n), work_tmp(2*n)

    ITYPE=1
    JOBZ='V'
    UPLO='U'
    LWORK=2*n

    call zhegv(ITYPE,JOBZ,UPLO,n,h,n,s,n,w,WORK,LWORK,RWORK,INFO)
    if (INFO /= 0) then
      write(*,*) 'ERROR: Generalized eigenvalue problem failed. zhegv failed with INFO =', INFO
      stop
    end if
  end

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         Bloch_Hamiltonian
!   INPUTS:       rkx,rky: Reciprocal space vector in fractional coordinates (with respect to the lattice vectors of the _tb.dat file)
!                 h_matrix: Coefficient matrices of the kp Hamiltonian polynomial expansion, read from Materialname_kp.dat (dimension n_tuples x norb x norb)
!                 norb: The dimension of the Hamiltonian matrix
!                 power_x,power_y: Powers of rkx and rky for each term of the kp Hamiltonian polynomial expansion
!   OUTPUTS:      hkernel: Bloch Hamiltonian matrix in reciprocal space
!                 skernel: Overlap matrix in reciprocal space
!                 akernel: Dipole matrix in reciprocal space
!                 hderkernel: Derivative of the Hamiltonian matrix
!                 sderkernel: Derivative of the overlap matrix
!   DESCRIPTION:  This subroutine is meant to evaluate kp Hamiltonians in reciprocal space. 
!                 It reads the kp Hamiltonian from the file Materialname_kp.dat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Bloch_Hamiltonian(rkx,rky,h_matrix,norb,power_x,power_y,hkernel,skernel,akernel,hderkernel,sderkernel)
    implicit none

    integer norb
    integer power_x(:),power_y(:)
    integer ii, i, n_tuples
    

    dimension skernel(norb,norb)
    dimension hkernel(norb,norb)
    dimension sderkernel(3,norb,norb)
    dimension hderkernel(3,norb,norb)
    dimension akernel(3,norb,norb)

    complex*16 skernel,sderkernel,hkernel,hderkernel,akernel
    complex*16 h_matrix(:,:,:)
    complex*16 :: Identity(norb,norb)


    real(8) rkx,rky,rkz

    hkernel=0.0d0
    hderkernel=0.0d0
    skernel=0.0d0
    sderkernel=0.0d0
    akernel=0.0d0

    Identity = (0.0d0, 0.0d0)
    do ii = 1, norb
      Identity(ii,ii) = (1.0d0, 0.0d0)
    end do
    
    n_tuples = size(power_x)

    hkernel = (0.0d0,0.0d0)

    do i = 1,n_tuples
      hkernel(:,:) = hkernel(:,:) + h_matrix(i,:,:)*rkx**power_x(i)*rky**power_y(i)
    end do

    hkernel(:,:) = hkernel(:,:) + Identity*(rkx**2 + rky**2)/2

    hderkernel = (0.0d0,0.0d0)
    do i = 1, n_tuples
      if (power_x(i) > 0) then
        hderkernel(1,:,:) = hderkernel(1,:,:) + power_x(i)*h_matrix(i,:,:)*rkx**(power_x(i)-1)*rky**power_y(i)
      end if
      if (power_y(i) > 0) then
        hderkernel(2,:,:) = hderkernel(2,:,:) + power_y(i)*h_matrix(i,:,:)*rkx**power_x(i)*rky**(power_y(i)-1)
      end if
    end do
    hderkernel(1,:,:) = hderkernel(1,:,:) + Identity*rkx
    hderkernel(2,:,:) = hderkernel(2,:,:) + Identity*rky
    hderkernel(3,:,:) = 0.0d0
    akernel = 0.0d0
    skernel(:,:)  = Identity(:,:)
    sderkernel = 0.0d0
    if (rkx**2 + rky**2 > 0.0005) then
      hkernel = (0.0d0,0.0d0)
      do ii = 1, norb
        hkernel(ii,ii) = dcmplx(dble(ii), 0.0d0)
      end do
      hderkernel = 0.0d0
    end if

  end subroutine 

end module constants_math

