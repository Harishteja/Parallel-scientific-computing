program main
implicit none
integer itr,n,i,j
double precision :: error,error_max,del,time_i,time_f
double precision ,allocatable :: x(:),y(:),phi(:,:),phi_old(:,:),q(:,:)

del =0.005
call cpu_time(time_i)
n=(2.0/del)+1
print*,'the n is',n
allocate(x(n),y(n))
allocate(phi(n,n),phi_old(n,n),q(n,n))

x(1)=-1.0

do i=2,n
    x(i)=x(i-1)+del
end do


y(1)=-1.0
do j=2,n
    y(j)=y(j-1)+del
end do

do j=1,n
    do i=1,n
        phi(i,j)=0.0
    end do
end do

phi_old=phi

do j=1,n
   phi(1,j)=sin(2.0*3.14*y(j))
end do

do i=1,n
phi(i,n)=0.0
phi(i,1)=0.0
enddo


do j=1,n
    do i=1,n
        q(i,j)=(x(i)**2)+(y(j)**2)
    end do
end do

error_max=0.0001
error=10

do while(error>error_max)
    error=0
    itr=itr+1
    do j=2,n-1
        do i=2,n-1
            phi(i,j)=(phi_old(i+1,j)+phi_old(i-1,j)+phi_old(i,j+1)+phi_old(i,j-1)+(del**2)*q(i,j))/4.0
        end do
    end do

    do j=1,n
    phi(n,j)=((4*phi(n-1,j))-phi(n-2,j))/3.0
    end do
do j=1,n
    do i=1,n
        error=error+abs(phi(i,j)-phi_old(i,j))
    end do
end do

do j=1,n
  do i=1,n
    phi_old(i,j)=phi(i,j)
  end do
end do 
print*,error
print*,itr
end do
call cpu_time(time_f)
print*,itr
print*,'the time taken is',time_f-time_i
OPEN(UNIT=100, FILE="serial_jacobi_x.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1,n
               WRITE(100,*) x(i)
         end do
OPEN(UNIT=100, FILE="serial_jacobi_y.dat", STATUS="REPLACE", ACTION="WRITE")

         do j = 1,n
               WRITE(100,*) y(j)
         end do
         
OPEN(UNIT=100, FILE="serial_jacobi_phi.dat", STATUS="REPLACE", ACTION="WRITE")

         do j = 1,n
            do i=1,n
               WRITE(100,*) phi(i,j)
            end do
         end do
OPEN(UNIT=100, FILE="phi_y.dat", STATUS="REPLACE", ACTION="WRITE")

         do j = 1,n
               WRITE(100,*) phi(101,j)
         end do
OPEN(UNIT=100, FILE="phi_x.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1,n
               WRITE(100,*) phi(i,101)
         end do

end program main




