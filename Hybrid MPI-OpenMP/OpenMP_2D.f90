program main
use omp_lib
    implicit none
    
    ! Input parameters
    integer ::itr, nx, ny, nt,i,j,threads
    real ::L,t1,t2, dt, alpha,time,dx,dy,k1,k2
    real, allocatable :: x(:),y(:)
    real, allocatable :: T(:,:), Told(:,:)
   
    threads=8
 !   allocate(T(nx,ny), Told(nx,ny))
!    allocate(x(nx),y(ny))

    t1=omp_get_wtime()
    L=256
    nx = 512 
    ny = 512 
    dt = 0.02 
    alpha = 1.4 
    nt = 2000 
    dx=L/nx
    dy=L/ny

    allocate(T(nx,ny), Told(nx,ny))
    allocate(x(nx),y(ny))
    
    x(1)=0
    y(1)=0
    
     do i=2,nx

       x(i)=x(i-1)+dx

     end do
     

     
     do i=2,ny

       y(i)=y(i-1)+dy

     end do
     
    
  T = 300
    do i = 1, nx
        T(1,i) = 600.0
        T(ny,i) = 900.0
    end do
   
    do j = 1, ny
        T(j,1) = 400.0
        T(j,nx) = 800.0
    end do
    
print*, 'okay'
    
!    T(1,1) = (T(1,2) + T(2,1)) / 2.0
 !   T(1,nx) = (T(1,nx-1) + T(2,nx)) / 2.0
 !   T(ny,1) = (T(ny-1,1) + T(ny,2)) / 2.0
  !  T(ny,nx) = (T(ny,nx-1) + T(ny-1,nx)) / 2.0
   

 
    Told=T
    k1 = alpha * (dt / (dx**2))
    k2 = alpha * (dt / (dy**2))

   time=0
   itr=0
do while (time< nt)
       itr=itr+1
       !$omp parallel num_threads(threads)

     !$omp do collapse(2) private(i,j) 
        do i = 2, nx-1
            do j = 2, ny-1
             
                    T(i,j)=Told(i,j) + (k1 * (Told(i+1,j) - 2.0*Told(i,j) + Told(i-1,j))) &
                                  + (k2 * (Told(i,j-1) - 2.0*Told(i,j) + Told(i,j+1)))

            end do
       end do
      !$omp end do

     !$omp end parallel 
      Told=T
        time=time+dt

 end do
t2=omp_get_wtime()
OPEN(UNIT=100, FILE="time_1.dat", STATUS="REPLACE", ACTION="WRITE")

 write(100,*) t2-t1
OPEN(UNIT=100, FILE="iterations_1.dat", STATUS="REPLACE", ACTION="WRITE")
write(100,*) itr
!do i=1,nx

 !  write(*,*) (T(i,j),j=1,ny)
   
!end do
!OPEN(UNIT=100, FILE="x.dat", STATUS="REPLACE", ACTION="WRITE")

 !        do i = 1,nx
!               WRITE(100,*) x(i)
  !       end do
!OPEN(UNIT=100, FILE="y.dat", STATUS="REPLACE", ACTION="WRITE")
  
  !         do i = 1,nx
 !                WRITE(100,*) y(i)
 !          end do


OPEN(UNIT=100, FILE="temperature_1.dat", STATUS="REPLACE", ACTION="WRITE")
 
             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,ny)
             end do

end program main


