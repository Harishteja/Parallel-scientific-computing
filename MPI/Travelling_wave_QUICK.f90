program main
implicit none
integer::i,x_min,n,n1,itr,j
real :: dt,dx,error,error_max,t,t_max
real,allocatable :: x(:),u(:),u_old(:),u_exact(:)
x_min=0

dx=0.002

n=int(2.0/dx)+1

allocate(x(n),u(n),u_old(n),u_exact(n))
x(1)=0.0
do i=2,n
  x(i)=x(i-1)+dx
end do
dt=0.0001
n1=int(0.5/dx)+1

do i=1,n
 if (x(i).le.0.5) then
   u(i)=sin(4*3.14*x(i))
 else 
   u(i)=0
 end if
end do
!do i = n1 + 1, n
!  u(i) = 0.0
!end do

!do j=1,n
!  if( (x(j)-1 ) <= 0.5) then
!            u_exact(j) = sin(4.0 * 3.14 * (x(j)-1 ))
!        else
!            u_exact(j) = 0.0 
!  end if
!end do
!do i=1,n
!u(i)=u_exact(i)
!end do
itr=0
t_max=1
t=0


do while(t<=t_max)
!t=t+dt
itr=itr+1
u(1)=0
u(n)=0
!u(2)=u_old(2)-(dt/dx)*(u_old(2)-u_old(1))
 do i=1,n
    u_old(i)=u(i)
 end do
!u(n)=0
 u(2)=u_old(2)-(dt/dx)*(u_old(2)-u_old(1))
 do i=3,n-1
   u(i)=u_old(i)-((dt/dx)*(((3.0/8.0)*u_old(i))-((7.0/8.0)*u_old(i-1))+((1.0/8.0)*u_old(i-2))+((3.0/8.0)*u_old(i+1))))
!     u(i)=u_old(i)-(dt/dx)*(u_old(i)-u_old(i-1))
end do

 do j=1,n
  if( ((x(j)-t) > 0.0).and.((x(j)-t ) <= 0.5)) then
            u_exact(j) = sin(4.0 * 3.14 * (x(j)-t ))
        else
            u_exact(j) = 0.0
  end if
end do
t=t+dt
end do


print *, "x        Numerical Solution        Exact Solution"
    do i = 1, n
        print *, x(i), u(i), u_exact(i)
    end do
    print *, "Number of iterations:", itr
    OPEN(UNIT=100, FILE="u_assgn_2.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1, n
               WRITE(100,*) u(i)
         end do
    OPEN(UNIT=100, FILE="u_exact_asign_2.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1,n
               WRITE(100,*) u_exact(i)
         end do
    OPEN(UNIT=100, FILE="x_i.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1, n
               WRITE(100,*) x(i)
         end do

end program main

