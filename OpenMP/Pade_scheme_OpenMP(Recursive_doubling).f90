program main
use omp_lib
implicit none
   integer:: i,j,k,thread
   real :: h,n,s,sumd,n_steps,m
   double precision:: t1,t2
   real,allocatable ::x(:),p(:),alpha(:),func(:),y(:),b(:),xf(:),e(:),g(:),f(:),beta(:),u(:),z(:)
   real, allocatable::a(:,:),lf(:,:),uf(:,:),check(:,:),check2(:),checkf(:)
   n=1000
   allocate(x(n),alpha(n),p(n),func(n),y(n),b(n),xf(n),g(n),e(n),f(n-1),z(n),check2(n),checkf(n))
   allocate(a(n,n),beta(n),u(n),lf(n,n),uf(n,n),check(n,n))
   h=3/(n-1)
   thread=8

   x(1)=0.00
   do i=2,n
     x(i)=x(i-1)+h
   end do

    t1=omp_get_wtime()
 !$omp parallel num_threads(thread)
   !$omp do
   do i=1,n
       func(i)=sin(5*x(i))
   end do
   !$omp end do

   !$omp do
   do i=1,n
       p(i)=5*cos(5*x(i))
  end do
  !$omp end do


 a(1,1)=1
 a(2,1)=1
 a(1,2)=2
 a(n,n)=1
 a(n,n-1)=2
 a(n-1,n)=1
 do i=2,n-1
   do j=2,n-1
     if(i.eq.j) then
       a(i,j)=4
       e(i)=a(i,j)
     else if(abs(i-j).eq.1) then
      a(i,j)=1
     else
     a(i,j)=0
     end if
   end do
 end do

 e(1)=A(1,1)
 e(n)=a(n,n)
 g(1)=0
 !$omp do
 do i=2,n-1
     g(i)=1
  end do
 !$omp end do

  g(n)=2
  f(1)=2
  !$omp do
  do i=2,n-1
     f(i)=1
  end do
  !$omp end do

  f(n)=0.0

 b(1)=((-2.5*func(1))+2*func(2)+0.5*func(3))/h
 !$omp do
 do i=2,n-1
   b(i)=(3*(func(i+1)-func(i-1)))/h
 end do
 !$omp end do
 b(n)=(2.5*func(n)-2*func(n-1)-0.5*func(n-2))/h

m=log(n)/log(2.0)
n_steps=ceiling(m)

do k=1,n_steps

 !$omp do private(i)
 do i=1,n
  if (i.ge.(2**(k-1) + 1)) then
     alpha(i)=-g(i)/e(i-2**(k-1))
     
   else 
    alpha(i)=0.0
    
    end if
  if(i.le.(n-2**(k-1))) then
   beta(i)=-f(i)/e(i+2**(k-1))
  else
   beta(i)=0.0
   end if
   if (i .ge. (2**k + 1)) then
       g(i) = alpha(i) * g(i-2**(k-1))
      else
       g(i) = 0.0
 end if
  if(i.le.(n-2**k)) then
   f(i)=beta(i)*f(i+2**(k-1))
   else
   f(i)=0.0
 end if
 e(i) = alpha(i)*f(i-2**(k-1))+e(i)+beta(i)*g(i+2**(k-1))
 b(i) = alpha(i)*b(i-2**(k-1))+b(i)+beta(i)*b(i+2**(k-1))

 end do
 !$omp end do
 end do


 !$omp do
 do i=1,n
   xf(i)=b(i)/e(i)
end do
!$omp end do


!$omp end parallel
t2=omp_get_wtime()

write(*,*)'thi finarl is',xf
write(*,*)'the time taken is',t2-t1
end program main
