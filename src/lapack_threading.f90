module lapack_threading
  ! OpenBLAS/MKL keep their own internal thread pool for LAPACK calls
  ! (zheev/zhegv). That pool is not safe to re-enter from several OpenMP
  ! threads at once: concurrent diagoz/diagoz_gen calls from different
  ! threads (in JDOS_sp, ome_sp, ...) can corrupt each other's LAPACK work.
  ! force_serial_lapack() pins the library to 1 internal thread, so every
  ! diagoz/diagoz_gen call stays serial and OpenMP k-point loops remain the
  ! only source of parallelism. This must only be active around those
  ! single-particle OMP-parallel diagonalization loops: call
  ! restore_parallel_lapack() right after them so LAPACK goes back to
  ! multithreaded for the rest of the code (exciton/BSE stages, ...).
  use iso_c_binding, only: c_int
  use omp_lib, only: omp_get_max_threads
  implicit none
  private
  public :: force_serial_lapack
  public :: restore_parallel_lapack

#ifdef USE_MKL
  interface
    subroutine blas_set_num_threads(num_threads) bind(C, name="MKL_Set_Num_Threads")
      import :: c_int
      integer(c_int), value :: num_threads
    end subroutine blas_set_num_threads
  end interface
#else
  interface
    subroutine blas_set_num_threads(num_threads) bind(C, name="openblas_set_num_threads")
      import :: c_int
      integer(c_int), value :: num_threads
    end subroutine blas_set_num_threads
  end interface
#endif

contains

  subroutine force_serial_lapack()
    call blas_set_num_threads(1)
  end subroutine force_serial_lapack

  subroutine restore_parallel_lapack()
    call blas_set_num_threads(omp_get_max_threads())
  end subroutine restore_parallel_lapack

end module lapack_threading
