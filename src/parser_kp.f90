module parser_kp
   use parser_input_file, only:read_line_numbers_int, iflag_orthonormal
   use parser_wannier90_tb, only: norb, R, material_name, nR, nRvec
   implicit none
   private
   public :: norb
   public :: R,h_matrix,power_x,power_y
   public :: material_name
   public :: kp_get

   integer i_tuple, max_order
   integer, allocatable :: power_x(:),power_y(:)

   complex*16, allocatable:: h_matrix(:,:,:)

contains

!--------------------------------------------------------------------
   pure function to_lower(str) result(out)
!> Convert a string to lower‑case (portable; no compiler extension)
      character(len=*), intent(in) :: str
      character(len=len(str))      :: out
      integer :: i
      out = str
      do i = 1, len(str)
         if (out(i:i) >= 'A' .and. out(i:i) <= 'Z') &
            out(i:i) = char(iachar(out(i:i)) + 32)
         end do
      end function to_lower
      !--------------------------------------------------------------------
      subroutine kp_get(material_name_in)
         implicit none
         ! -----------------------------------------------------------------
         character(len=*), intent(in) :: material_name_in    ! full path read from input.txt
         ! -----------------------------------------------------------------
         integer          :: fp, i, ialpha, ialphap, n_tuples
         integer          :: nkk1, nkk2
         real(8)          :: a1,a2,a3,a4,a5,a6
         character(len=:), allocatable :: file2open, basename
         integer          :: p, ext_pos
         integer :: num_chunks
         integer :: rem
         ! -----------------------------------------------------------------
         write(*,*) '2. Entering parser_kp'

      ! === 1.  Use the path exactly as supplied ========================
      !write(*,*) "MATERIAL NAME PARSED:", material_name_in
      file2open = trim(material_name_in)

      ! === 2.  Derive clean material name (no dir, no _tb.dat) =========
      p = max( index(file2open,'/',back=.true.),  &
         index(file2open,'\',back=.true.) )   ! works on Win/Linux
      if (p == 0) then
         basename = file2open
      else
         basename = file2open(p+1:)
      end if

      ext_pos = len_trim(basename) - len('_kp.dat') + 1
      if ( ext_pos > 0 .and. to_lower(basename(ext_pos:)) == '_kp.dat' ) then
         basename = basename(:ext_pos-1)
      end if
      material_name = adjustl(basename)

      ! === 3.  Open kp Hamiltonian file ==================================
      open(unit=fp, file=file2open, action='read', status='old')
      read(fp,*)
      read(fp,*) R(1,1),R(1,2),R(1,3)
      read(fp,*) R(2,1),R(2,2),R(2,3)
      read(fp,*) R(3,1),R(3,2),R(3,3)
      read(fp,*) norb
      read(fp,*) max_order

      n_tuples = (max_order+1)*(max_order+2)/2
      allocate(power_x(n_tuples), power_y(n_tuples))
      allocate(h_matrix(n_tuples, norb, norb))

      !get the Hamiltonian matrix
      do i_tuple =1,(max_order+1)*(max_order+2)/2  ! Number of permutations (n,m) where n,m are kx^nk_y^m for a maximum order.
         read(fp,*) power_x(i_tuple),power_y(i_tuple)
         do ialphap=1,norb
            do ialpha=1,norb
               read(fp,*) nkk1,nkk2,a1,a2
               h_matrix(i_tuple,nkk1,nkk2)=complex(a1,a2)
            end do
         end do
         read(fp,*)
      end do

      close(fp)

      !convert units: to Hartree and bohrs
      !hhop=hhop/27.211385d0
      !rhop_c=rhop_c/0.52917721067121d0
      R=R/0.52917721067121d0
      nR = 0
      ! The following makes the system 2D in x and y
      allocate(nRvec(2,3))
      nRvec = 0
      nRvec(1,1) = 1
      nRvec(2,2) = 1
      write(*,*) '   kp Hamiltonian has been read'
   end subroutine kp_get


end module parser_kp
