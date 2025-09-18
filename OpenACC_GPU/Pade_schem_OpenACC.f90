program main
use openacc
implicit none
   integer:: i,j,k
   real :: h,s,sumd,n
   double precision :: start,end_time
  real,allocatable ::x(:),p(:),func(:),y(:),b(:),xf(:),e(:),g(:),f(:),l(:),u(:),z(:)
   real, allocatable::a(:,:),lf(:,:),uf(:,:)
   n=1000
  
   allocate(x(n),p(n),func(n),y(n),b(n),xf(n),g(n-1),e(n),f(n-1),z(n))
   allocate(a(n,n),l(n-1),u(n),lf(n,n),uf(n,n))
   h=3/(n-1)
call cpu_time(start)
   x(1)=0.00
   do i=2,n
     x(i)=x(i-1)+h
   end do
 !$acc  parallel loop num_gangs(100)
   do i=1,n
       func(i)=sin(5*x(i))
       p(i)=5*cos(5*x(i))
  end do
!$acc end parallel loop


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

 do i=1,n-2
     g(i)=1
  end do

  g(n-1)=2
  f(1)=2
  do i=2,n-1
     f(i)=1
  end do

 b(1)=((-2.5*func(1))+2*func(2)+0.5*func(3))/h

 do i=2,n-1
   b(i)=(3*(func(i+1)-func(i-1)))/h
 end do

 b(n)=(2.5*func(n)-2*func(n-1)-0.5*func(n-2))/h
 
 u(1)=e(1)
 do i=1,n-1
    l(i)=g(i)/u(i)
    u(i+1)=e(i+1)-(l(i)*f(i))
 end do

 do i=1,n
   do j=1,n
    if (i.eq.j) then
     lf(i,j)=1
    else if((i-j).eq.1) then
      lf(i,j)=l(j)
     else
     lf(i,j)=0
     end if
   end do
  end do

  do i=1,n
     do j=1,n
      if (i.eq.j) then
       uf(i,j)=u(j)
      else if((j-i).eq.1) then
        uf(i,j)=f(i)
       else
       uf(i,j)=0
       end if
    end do
  end do

 z(1)=b(1)/lf(1,1)
!$acc  loop seq
do i=2,n
    sumd=0
 do j=1,i-1
    sumd=sumd+(lf(i,j)*z(j))
 end do
 z(i)=b(i)-sumd
end do
!$acc end  loop

xf(n)=z(n)/uf(n,n)
!$acc  loop seq
do i=n-1,1,-1
 sumd=0
 do j=i+1,n
    sumd=sumd+(uf(i,j)*xf(j))
 end do
 xf(i)=(z(i)-sumd)/uf(i,i)
end do
!$acc end loop

call cpu_time(end_time)
write(*,*)'the final answer x is',xf
write(*,*)'time taken is',end_time-start
OPEN(UNIT=100, FILE="x.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1,n
               WRITE(100,*) xf(i)
         end do
end program main
 

