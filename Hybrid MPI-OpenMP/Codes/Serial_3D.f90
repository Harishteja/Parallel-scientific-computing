program main
    implicit none

    ! Input parameters
    integer ::nz, nx, ny, nt,i,j,k,itr
    real ::L,t1,t2, dt, alpha,time,dx,dy,dz,k1,k2,k3
    real, allocatable :: x(:),y(:),z(:)
    real, dimension(:,:,:), allocatable :: T, Told

   call cpu_time(t1)
   L=100    
    nx =64 
    ny =64
    nz =64
    dt = 0.02 
    alpha = 1.4 
    nt = 2000 
    dx=L/nx
    dy=L/ny
    dz=L/nz
 
    allocate(T(nx,ny,nz), Told(nx,ny,nz))
    allocate(x(nx),y(ny),z(nz))
    x(1)=0
    y(1)=0
    z(1)=0
     do i=2,nx

       x(i)=x(i-1)+dx

     end do
     do i=2,ny

       y(i)=y(i-1)+dy

     end do
     do i=2,nz
        
       z(i)=z(i-1)+dz

     end do
    
    T = 300.0
    
    do i = 1, nx
        T(1,i,:) = 600.0
        T(ny,i,:) = 900.0
    end do
    do j = 1, ny
        T(j,1,:) = 400.0
        T(j,nx,:) = 800.0
    end do
    T(:,:,1)=700
    T(:,:,nz)=500
!    T(1,1) = (T(1,2) + T(2,1)) / 2.0
!    T(1,nx) = (T(1,nx-1) + T(2,nx)) / 2.0
!    T(ny,1) = (T(ny-1,1) + T(ny,2)) / 2.0
!    T(ny,nx) = (T(ny,nx-1) + T(ny-1,nx)) / 2.0

   
   ! x = linspace(0.0, 1.0, nx)
   ! y = linspace(0.0, 1.0, ny)
   

    Told = T

    ! CFL values for stability
   ! real :: k1, k2, CFL
    k1 = alpha * (dt / (dx**2))
    k2 = alpha * (dt / (dy**2))
    k3 = alpha * (dt / (dz**2))
   ! CFL = k1 + k2
   time=0
   itr=0
    ! Time loop
do while (time<nt )
itr=itr+1
        do i = 2, nx-1
            do j = 2, ny-1
               do k = 2, nz-1
                T(i,j,k) = Told(i,j,k) + (k1 * (Told(i+1,j,k) - 2.0*Told(i,j,k) + Told(i-1,j,k))) &
                                  + (k2 * (Told(i,j-1,k) - 2.0*Told(i,j,k) + Told(i,j+1,k))) &
                                  + k3 * (Told(i,j,k+1) - 2.0*Told(i,j,k) + Told(i,j,k-1))
            end do
        end do
      end do
        Told = T
        time=time+dt
 end do
call cpu_time(t2)
OPEN(UNIT=100, FILE="time_3d_serial.dat", STATUS="REPLACE", ACTION="WRITE")
 write(100,*) t2-t1
OPEN(UNIT=100, FILE="iterations_3d.dat", STATUS="REPLACE", ACTION="WRITE")
 write(100,*) itr
 print*, 'the time is',time
do i=1,nx

   write(*,*) (T(i,2,j),j=1,ny)
   
end do
!OPEN(UNIT=100, FILE="x.dat", STATUS="REPLACE", ACTION="WRITE")

 !        do i = 1,nx
   !            WRITE(100,*) x(i)
  !       end do
!OPEN(UNIT=100, FILE="y.dat", STATUS="REPLACE", ACTION="WRITE")
  
  !         do i = 1,nx
   !              WRITE(100,*) y(i)
    !       end do


OPEN(UNIT=100, FILE="temperature_3d_1.dat", STATUS="REPLACE", ACTION="WRITE")
 
             do i = 1,nx
                 WRITE(100,*) (T(i,32,k),k=1,ny)
             end do

end program main


