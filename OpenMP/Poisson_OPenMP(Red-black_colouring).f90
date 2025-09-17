program main
use omp_lib
implicit none
    integer :: i,n, j,p,threads,nd,e,istart,iend
    real :: m
    double precision :: error,errormax,del,l,c,d,error3
    double precision, dimension(:),allocatable :: x,y
    double precision, dimension(:,:),allocatable ::phi,phi_old,phi_a,q,error2,k
    del=0.005
    m=2.0/del
    n=m+1
    threads=16
    write(*,*),n
    ALLOCATE(y(n))
    allocate(x(n))
    allocate(phi(n,n),q(n,n),error2(n,n))
    allocate(phi_old(n,n),phi_a(n,n),k(n-2,n-2))
    
   !$omp parallel num_threads(threads)
    c=omp_get_wtime()
    
     x(1)=-1
    do i=2,n
       x(i)=x(i-1)+del
    end do
    
    y(1)=-1
    do j=2,n
       y(j)=y(j-1)+del
    end do
    
    !$omp do collapse(2)
    do i=1,n
      do j=1,n
        phi(i,j)=0
      end do
   end do
   !$omp end do

   !$omp do collapse(2)
   do i=1,n
      do j=1,n
        q(i,j)=2*(2-(x(i)*x(i))-(y(j)*y(j)))
      end do
   end do 
   !$omp end do
 
   errormax=0.1
   error=10
   p=0

   !$omp do collapse(2)
   do i=1,n
     do j=1,n
        phi_a(i,j)=(x(i)*x(i)-1)*(y(j)*y(j)-1)
      end do
   end do
   !$omp end do

 !$omp end parallel

  do while (error>errormax)
     p=p+1
     error=0
     nd=2*n-1
 
 !$omp parallel num_threads(threads)

      do e=1,nd
    
          if(e.le.n) then
              istart=1
              iend=e
          else
             istart=e-n+1
             iend=n
          end if
       
     !$omp do private(i,j)
      do i=istart,iend
         j=e-i+1
         if ((i .gt. 1) .and. (i .lt. N) .and. (j .gt. 1) .and. (j .lt. N)) then
         phi(i,j)=(phi(i+1,j)+phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+(del*del*q(i,j)))/4
         end if
      end do
   !$omp end do
     end do
!$omp end parallel
     do i=2,n-1
        do j=2,n-1
           error =error+ABS((phi_a(i,j)-phi(i,j))/phi_a(i,j))
        end do
     end do
end do
 d= omp_get_wtime()
      print*,'time taken',d-c
      write(*,*)'no of iterations',p
 end program main
