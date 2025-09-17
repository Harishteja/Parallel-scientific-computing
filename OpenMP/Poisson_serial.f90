program main
implicit none
    integer :: i, j, n,p
    real :: m,t2,t1
    double precision :: error,errormax,del
    double precision, dimension(:),allocatable :: x,y
    double precision, dimension(:,:),allocatable ::phi,phi_old,phi_a,q
    del=0.1
    m=2.0/del
    n=m+1
    ALLOCATE(y(n))
    allocate(x(n))
    allocate(phi(n,n),q(n,n))
    allocate(phi_old(n,n),phi_a(n,n))
 
    call cpu_time(t1)

    x(1)=-1
    do i=2,n
       x(i)=x(i-1)+del
    end do

    y(1)=-1
    do j=2,n
       y(j)=y(j-1)+del
    end do

    do i=1,n
      do j=1,n
        phi(i,j)=0
      end do
    end do
 
   do i=1,n
      do j=1,n
        q(i,j)=2*(2-(x(i)*x(i))-(y(j)*y(j)))
      end do
   end do  
   
   errormax=0.01
   error=10
   p=0

   do i=1,n
     do j=1,n
        phi_a(i,j)=(x(i)*x(i)-1)*(y(j)*y(j)-1)
      end do
   end do

do while (error>errormax)
     p=p+1
     error = 0
     phi_old=phi
     do j=2,n-1
       do i=2,n-1
          phi(i,j)=(phi_old(i+1,j)+phi(i-1,j)+phi_old(i,j+1)+phi(i,j-1)+(del*del*q(i,j)))/4
       end do
     end do

     do i=2,n-1
        do j=2,n-1
           error = error+ABS((phi_a(i,j)-phi(i,j))/phi_a(i,j))
        end do
    end do
 end do
  call cpu_time(t2)

      write(*,*)'no of iterations',p
      write(*,*)'time taken is',t2-t1
 end program main
