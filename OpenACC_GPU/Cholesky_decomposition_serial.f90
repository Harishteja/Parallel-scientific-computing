program cholesky_decomposition
  implicit none
  

  integer, parameter :: N = 10
  real, dimension(N,N) :: a, l
  double precision :: start_time,end_time
  real :: cpu_time_used
  integer :: i, j, k
 call cpu_time(start_time)
  call init_matrix(a)
  l = 0.0
  call cholesky(a, l)
  call print_matrix(l)
  
  call cpu_time(cpu_time_used)
  call cpu_time(end_time)
  print *, "Time taken by the serial code for N =", N, "is :", end_time-start_time, "seconds"
        
contains

 
  subroutine init_matrix(mat)
    real, dimension(N,N) :: mat
    integer :: i, j
    
    do i = 1, N
      do j = 1, i
        mat(i,j) = (i + j) / (N * 1.0) / N
        mat(j,i) = mat(i,j)
      end do
      mat(i,i) = 1.0
    end do
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
    
    do i = 1, N
      do j = 1, i
        sum = 0.0
        do k = 1, j-1
          sum = sum + l(i,k) * l(j,k)
        end do
        if (i == j) then
          l(i,j) = sqrt(a(i,i) - sum)
        else
          l(i,j) = (a(i,j) - sum) / l(j,j)
        end if
      end do
    end do
  end subroutine cholesky


end program cholesky_decomposition


