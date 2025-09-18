program cholesky_decomposition
use openacc
  implicit none
  integer, parameter :: N = 10
  real, dimension(N,N) :: a, l
  double precision :: start_time,end_time
  real :: cpu_time_used
  integer :: i, j, k

    ! do i=1,n
   !  do j=1,n
  !     l(i,j)=0
 !   end do
!    end do
call cpu_time(start_time)
!!$acc enter data create(a(0:N)(0:N)), copy(L[0:N][0:N]) 
!$acc  data create(a(1:N,1:N)), copy(L(1:N,1:N))

  call init_matrix(a)
  l=0.0
  call cholesky(a, l)
   
  ! !call print_matrix(l)
  ! call cpu_time(cpu_time_used)
  ! print *, "Time taken by the serial code for N =", N, "is :", cpu_time_used, "seconds"
  ! !$acc exit data
  
!$acc end data
!call print_matrix(L)
call cpu_time(end_time)
 call cpu_time(cpu_time_used)
  print *, "Time taken by the serial code for N =", N, "is :", end_time-start_time, "seconds"
contains

  subroutine init_matrix(mat)
    real, dimension(N,N) :: mat
    integer :: i, j
    !$acc parallel loop present(mat)
    do i = 1, N
      do j = 1, i-1
        mat(i,j) = (i + j) / (N * 1.0) / N
        mat(j,i) = mat(i,j)
      end do
      mat(i,i) = 1.0
    end do
    !$acc end parallel
  end subroutine init_matrix
  
  subroutine print_matrix(mat)
     real, dimension(N,N) :: mat
     integer :: i, j

     do i = 1, N
       do j = 1, N
         write(*, '(f8.2)', advance='no') mat(i,j)
       end do
       write(*, *)
     end do
   end subroutine print_matrix

 
  subroutine cholesky(a, l)
    real, dimension(N,N) :: a, l
    integer :: i, j, k
    real :: sum
!!$acc parallel loop present(a) num_gangs(1) 
   do i = 1, N
  !$acc parallel loop present(a(1:N,1:N),L(1:N,1:N)) num_gangs(1)
  do j = 1, i
        sum = 0.0
        !$acc loop reduction(+:sum)
        do k = 1, j-1
          sum = sum + l(i,k) * l(j,k)
        end do
        !$acc end loop
        if (i == j) then
          l(i,j) = sqrt(a(i,i) - sum)
        else
          l(i,j) = (a(i,j) - sum) / l(j,j)
        end if
      end do
      !$acc end parallel loop
    end do
   ! !$acc end parallel loop
 end subroutine cholesky


end program cholesky_decomposition

