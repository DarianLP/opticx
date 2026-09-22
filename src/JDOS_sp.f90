module JDOS_sp
  use constants_math
  use lapack_threading, only:force_serial_lapack,restore_parallel_lapack
  use parser_input_file, &
  only:e1,e2,eta,nw,broadening_type_text,iflag_orthonormal
  use parser_wannier90_tb, &
  only:material_name,norb,nR, nRvec, R, shop, hhop
  use parser_optics_xatu_dim, &
  only:npointstotal,vcell, &
  nv_ex,nband_ex,nband_index, &
  rkxvector,rkyvector,rkzvector !k-vectors only used for testing
  use omp_lib

  implicit none

  contains

  subroutine get_JDOS_sp()
    implicit none
    integer :: ibz
    integer :: i,ii

    ! Per-thread scratch arrays (listed in PRIVATE clause below)
    real(8)     :: e_nband_local(nband_ex)
    real(8)     :: rkx,rky,rkz
    real(8)     :: e(norb)
    complex(8)  :: hkernel(norb,norb)
    complex(8)  :: skernel(norb,norb)

    ! frequency grid and conductivity tensor
    real(8)     :: wp(nw)
    real(8)     :: eta1
    real(8) :: JDOS_sp(nw)

    write(*,*) 'Computing JDOS_sp'

    !initialize sp arrays
    JDOS_sp = 0.0d0
    !build frequency grid (Hartree) and broadening (Hartree)
    call fill_frequency_grid(nw,e1,e2,eta,wp,eta1)
    ! ----------------------------------------------------------------
    !
    ! * PRIVATE  : each thread owns its own scratch arrays.
    ! * REDUCTION: JDOS_w_sp is accumulated across threads safely.
    ! ----------------------------------------------------------------
    ! diagoz/diagoz_gen are called concurrently below, so LAPACK must stay
    ! pinned to 1 thread for the duration of this loop (see
    ! lapack_threading); restored right after.
    call force_serial_lapack()
    !$OMP PARALLEL DO          &
    !$OMP   DEFAULT(SHARED)    &
    !$OMP   PRIVATE(ibz, i, ii, rkx, rky, rkz, hkernel, skernel, e, e_nband_local) &
    !$OMP   REDUCTION(+:JDOS_sp) &
    !$OMP   SCHEDULE(dynamic)
    do ibz=1,npointstotal
      rkx=rkxvector(ibz)
      rky=rkyvector(ibz)
      rkz=rkzvector(ibz)
      call get_hk_sk(rkx,rky,rkz,norb,hkernel,skernel)
      e=0.0d0
      if (iflag_orthonormal) then
        call diagoz(norb,e,hkernel)
      else
        call diagoz_gen(norb,e,hkernel,skernel)
      end if

      !keep only the bands requested in the bandlist
      do i=1,nband_ex
        ii=nband_index(i)
        e_nband_local(i)=e(ii)
      end do

      !fill JDOS(w) for a given k point
      call evaluate_JDOS_sp(nband_ex,npointstotal,vcell,e_nband_local,nw,wp,eta1,JDOS_sp)
    end do
    !$OMP END PARALLEL DO
    call restore_parallel_lapack()

    !print JDOS
    write(*,*) '   Printing JDOS...'
    call print_JDOS_sp(nw,wp,JDOS_sp)

  end subroutine get_JDOS_sp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_hk_sk(rkx,rky,rkz,norb,hkernel,skernel)
    implicit none
    integer, intent(in) :: norb
    integer :: ialpha,ialphap,iRp
    real(8), intent(in) :: rkx,rky,rkz
    real(8) :: Rx,Ry,Rz
    complex(8), intent(out) :: hkernel(norb,norb),skernel(norb,norb)
    complex(8) :: phase,factor

    hkernel = 0.0d0
    skernel = 0.0d0

    do ialpha=1,norb
      do ialphap=1,ialpha
        do iRp=1,nR
          Rx=dble(nRvec(iRp,1))*R(1,1)+dble(nRvec(iRp,2))*R(2,1)+dble(nRvec(iRp,3))*R(3,1)
          Ry=dble(nRvec(iRp,1))*R(1,2)+dble(nRvec(iRp,2))*R(2,2)+dble(nRvec(iRp,3))*R(3,2)
          Rz=dble(nRvec(iRp,1))*R(1,3)+dble(nRvec(iRp,2))*R(2,3)+dble(nRvec(iRp,3))*R(3,3)
          phase=complex(0.0d0,rkx*Rx+rky*Ry+rkz*Rz)
          factor=exp(phase)
          hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+ &
          factor*hhop(iRp,ialpha,ialphap)
          skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
          factor*shop(iRp,ialpha,ialphap)
        end do
        hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
        skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))
      end do
    end do
  end subroutine get_hk_sk
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine print_JDOS_sp(nw,wp,JDOS_sp)
    implicit none
    integer,    intent(in) :: nw
    integer                :: iw
    real(8),    intent(in) :: wp(nw)
    real(8),intent(in) :: JDOS_sp(nw)
    real(8)                :: feps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    !write frequency dependent JDOS	  
    open(50,file='JDOS_sp_'//trim(material_name)//'.dat')

    do iw=1,nw
      feps=1.0d0 !use atomic units
      write(50,*) wp(iw)*27.211385d0, &
        feps*JDOS_sp(iw)
    end do
    close(50)

  end subroutine print_JDOS_sp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine evaluate_JDOS_sp(nband_ex,npointstotal,vcell,e,nw,wp,eta1,JDOS_sp)
    implicit none
    integer,    intent(in)    :: nw, nband_ex, npointstotal
    integer                   :: iw, nn, nnp
    real(8),    intent(in)    :: e(nband_ex), wp(nw), eta1, vcell
    real(8), intent(inout) :: JDOS_sp(nw)

    real(8)    :: fnn, fnnp, factor1
    real(8)    :: delta_nnp
    real(8)    :: JDOS_local(nw)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    JDOS_local = 0.0d0

    do nn=1,nband_ex
      !fermi distribution
      if (nn.le.nv_ex) then
        fnn=1.0d0
      else
        fnn=0.0d0
      end if

      do nnp=1,nband_ex
        !fermi distribution
        if (nnp.le.nv_ex) then
          fnnp=1.0d0
        else
          fnnp=0.0d0
        end if

        ! added check for degeneracy points in energy to prevent division by zero in the denominators
        if (abs(fnn-fnnp).lt.0.1d0 .or. abs(e(nn)-e(nnp)).lt.1.0d-6) then
          factor1=0.0d0
        else
          factor1=fnn-fnnp
        end if

        ! Skip entirely if no contribution
        if (factor1 == 0.0d0) cycle

        do iw=1,nw
          if (trim(broadening_type_text) == 'gaussian') then
            delta_nnp = -1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta1**2)*(wp(iw)-e(nn)+e(nnp))**2)
          else if (trim(broadening_type_text) == 'lorentzian') then
            delta_nnp = 1.0d0/pi*aimag(1.0d0/(-wp(iw)+e(nn)-e(nnp)+&
              complex(0.0d0,eta1)))
          else
            delta_nnp = -1.0d0/eta1*1.0d0/sqrt(2.0d0*pi)*&
              exp(-0.5d0/(eta1**2)*(wp(iw)-e(nn)+e(nnp))**2)
          end if

        JDOS_local(iw) = JDOS_local(iw) + &
                1/(dble(npointstotal)*vcell)*factor1*delta_nnp

        end do ! iw
      end do   ! nnp
    end do     ! nn

    ! Accumulate into the caller's array (protected by the outer
    ! k-point REDUCTION in get_JDOS_sp)
    JDOS_sp = JDOS_sp + JDOS_local

  end subroutine evaluate_JDOS_sp
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fill_frequency_grid(nw,e1,e2,eta,wp,eta1)
    implicit none
    integer, intent(in)  :: nw
    integer              :: i
    real(8), intent(in)  :: e1,e2,eta
    real(8), intent(out) :: wp(nw), eta1
    real(8)              :: wrange

    wp=0.0d0
    wrange=e2-e1
    do i=1,nw
      wp(i) = (e1+wrange/dble(nw)*dble(i-1))/27.211385d0
    end do
    eta1 = eta/27.211385d0
  end subroutine fill_frequency_grid
end module JDOS_sp