program main
use mpi
use omp_lib
implicit none

    ! Input parameters
    integer ::itr,threads, nx, ny, nt,i,j,my_rank,num_procs,p,start_idx,end_idx,ierr,v
    real ::L,dt, alpha,time,dx,dy,k1,k2,y_min
    double precision :: t1,t2
    real, allocatable :: x(:),y(:),T_c(:),T_d(:)
    real, dimension(:,:), allocatable :: T, Told,all_T

  
   ! nx = 12
   ! ny = 12 
  !  dt = 0.002 
 !   alpha = 1.4 
!    nt = 500 
!    dx=1
!    dy=1

!    allocate(T(nx,ny), Told(nx,ny))
  !  allocate(x(nx),y(ny))
call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierr)
t1=mpi_wtime()
p=num_procs

threads=8
L=256
nx = 512
    ny = 512
    dt = 0.02
    alpha = 1.4
    nt = 2000
    dx=L/nx
    dy=L/ny
start_idx = 1 + (my_rank * ny) / p
end_idx = (my_rank + 1) * ny / p

allocate(y(end_idx - start_idx + 1),x(nx),T_c(nx),T_d(nx))
   x(1)=0
   y_min=0
!    y(1)=0
     do i=2,nx

       x(i)=x(i-1)+dx

     end do

     do i = start_idx, end_idx
          y(i - start_idx + 1) = y_min + (i - 1) * dy
     end do

!     do i=2,ny

 !      y(i)=y(i-1)+dy

  !   end do
     if(my_rank==0) then
        print*,'harish',start_idx
        print*,'harish',end_idx

     end if
    allocate(T(nx,end_idx-start_idx+1))
    allocate(Told(nx,end_idx-start_idx+1))
    allocate(all_T(nx, ny))
   print*,'harish22'

   T = 300.0
    
    do i = start_idx,end_idx
        T(1,i-start_idx+1) = 600.0
        T(ny,i-start_idx+1) = 900.0
    end do
if (my_rank==0) then
    do i = 1, nx
        T(i,1) = 400.0
      !  T(j,nx) = 800.0
    end do
end if

if (my_rank==num_procs-1) then
  do i=1,nx
      T(i,end_idx-start_idx+1)=800.0
  end do
end if

if (my_rank==1) then
        do i = 1,nx
 !           WRITE(*,*) 'the temparaturein rannk 0 is',(T(i,j),j=1,end_idx-start_idx+1)
         end do
    end if
    Told=T
     if (my_rank==2) then
        do i = 1,nx
  !          WRITE(*,*) 'the temparature old in rannk 0 is',(Told(i,j),j=1,end_idx-start_idx+1)
         end do
    end if

   ! T(1,1) = (T(1,2) + T(2,1)) / 2.0
   ! T(1,nx) = (T(1,nx-1) + T(2,nx)) / 2.0
   ! T(ny,1) = (T(ny-1,1) + T(ny,2)) / 2.0
   ! T(ny,nx) = (T(ny,nx-1) + T(ny-1,nx)) / 2.0

   
   ! x = linspace(0.0, 1.0, nx)
   ! y = linspace(0.0, 1.0, ny)
   

 !   Told = T

    ! CFL values for stability
   ! real :: k1, k2, CFL
    k1 = alpha * (dt / (dx**2))
    k2 = alpha * (dt / (dy**2))
   ! CFL = k1 + k2
   time=0;
   itr=0;
    ! Time loop
do while( time<nt)
     if (my_rank < num_procs-1) then
        call MPI_Send(T(1:nx,end_idx-start_idx+1), nx, mpi_real, my_rank+1, 0, mpi_comm_world, ierr)
        call MPI_Recv(T_d, nx, mpi_real, my_rank+1, 0, mpi_comm_world, MPI_STATUS_IGNORE, ierr)

     end if
     if(my_rank>0) then

           call mpi_Send(T(1:nx,1), nx, mpi_real, my_rank-1, 0, mpi_comm_world, ierr)
           call mpi_Recv(T_c, nx, mpi_real, my_rank-1, 0, mpi_comm_world, MPI_STATUS_IGNORE, ierr)

     end if
    !$omp parallel num_threads(threads) private(i,j)
    v=omp_get_thread_num()
!     print*,'hello from harishteja from thread',v,'from processor',my_rank
       if(my_rank==0) then

         !$omp do collapse(2)
          do i = 2, nx-1
            do j = start_idx+1, end_idx-1
                T(i,j-start_idx+1) = Told(i,j-start_idx+1) + (k1 * (Told(i+1,j-start_idx+1) - 2.0*Told(i,j-start_idx+1)&
                                    + Told(i-1,j-start_idx+1))) &
                                  + (k2 * (Told(i,j-start_idx+1-1) - 2.0*Told(i,j-start_idx+1) + Told(i,j-start_idx+1+1)))
            end do
        end do
        !$omp end do

        !$omp do
         do i = 2, nx-1

                T(i,end_idx-start_idx+1) = Told(i,end_idx-start_idx+1) + (k1 * (Told(i+1,end_idx-start_idx+1)&
                                   - 2.0*Told(i,end_idx-start_idx+1) + Told(i-1,end_idx-start_idx+1))) &
                                  + (k2 * (Told(i,end_idx-start_idx+1-1) - 2.0*Told(i,end_idx-start_idx+1) + T_d(i)))
            end do
	    !$omp end do

     
     else if (my_rank==num_procs-1) then

      !$omp do collapse(2)
          do i = 2, nx-1
            do j = start_idx+1,end_idx-1
                T(i,j-start_idx+1) = Told(i,j-start_idx+1) + (k1 * (Told(i+1,j-start_idx+1) - 2.0*Told(i,j-start_idx+1)&
                                  + Told(i-1,j-start_idx+1))) &
                                  + (k2 * (Told(i,j-start_idx+1-1) - 2.0*Told(i,j-start_idx+1) + Told(i,j-start_idx+1+1)))
            end do
        end do
       !$omp end do
       !$omp do
         do i=2,nx-1
                
               T(i,1)=Told(i,1) + (k1 * (Told(i+1,1) - 2.0*Told(i,1) + Told(i-1,1))) &
                                  + (k2 * (T_c(i) - 2.0*Told(i,1) + Told(i,1+1)))
         end do
       !$omp end do

else

        !$omp do collapse(2)
         do i = 2, nx-1
            do j=start_idx+1,end_idx-1
                T(i,j-start_idx+1) =  Told(i,j-start_idx+1) + (k1 * (Told(i+1,j-start_idx+1)&
                                    - 2.0*Told(i,j-start_idx+1) + Told(i-1,j-start_idx+1))) &
                                    + (k2 * (Told(i,j-start_idx+1-1)- 2.0*Told(i,j-start_idx+1) + Told(i,j-start_idx+1+1)))
            end do
        end do
       !$omp end do

       !$omp do
         do i=2,nx-1
            
             T(i,1)=Told(i,1) + (k1 * (Told(i+1,1) - 2.0*Told(i,1) + Told(i-1,1))) &
                             + (k2 * (T_c(i) - 2.0*Told(i,1) + Told(i,1+1)))
         end do
         !$omp end do

	 !$omp do
        do i=2,nx-1
            
              T(i,end_idx-start_idx+1) = Told(i,end_idx-start_idx+1) + (k1 * (Told(i+1,end_idx-start_idx+1)&
                                       - 2.0*Told(i,end_idx-start_idx+1) + Told(i-1,end_idx-start_idx+1))) &
                                       + (k2 * (Told(i,end_idx-start_idx+1-1) - 2.0*Told(i,end_idx-start_idx+1) + T_d(i)))
         end do
	 !$omp end do

    end if
       !$omp end parallel
        Told = T
        time=time+dt
     
if (my_rank==1) then
 print*,'iteration is',itr
end if

end do
t2=mpi_wtime()
print*,'the time taken is',t2-t1

  if (my_rank==0) then
        do i = 1,nx
!            WRITE(*,*) 'the temparature final in rank 2 is',(T(i,j),j=1,end_idx-start_idx+1)
         end do
    end if

     call MPI_Allgather(T, nx*ny/p, mpi_real, all_T, nx*ny/p, mpi_real, MPI_COMM_WORLD, ierr)

 print*, 'the time is',time
 if (my_rank==0) then
 do i=1,nx
!   write(*,*) 'the temp is',(all_T(i,j),j=1,ny )
   
end do
end if
print*,'itr is',itr
OPEN(UNIT=100, FILE="time_hybrid.dat", STATUS="REPLACE", ACTION="WRITE")
   WRITE(100,*) t2-t1
OPEN(UNIT=100, FILE="procs_hybrid.dat", STATUS="REPLACE", ACTION="WRITE")
   WRITE(100,*) p
!         do i = 1,nx
 !              WRITE(100,*) x(i)
  !       end do
!OPEN(UNIT=100, FILE="y.dat", STATUS="REPLACE", ACTION="WRITE")
  
   !        do i = 1,nx
    !             WRITE(100,*) y(i)
     !      end do
OPEN(UNIT=100, FILE="temperature_hybrid_2.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (all_T(i,j),j=1,ny)
             end do



if (my_rank==0) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_0.dat", STATUS="REPLACE", ACTION="WRITE")
 
             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==1) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_1.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==2) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_2.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==3) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_3.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==4) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_4.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==5) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_5.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==6) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_6.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if
if (my_rank==7) then
OPEN(UNIT=100, FILE="temperature_hybrid_2_7.dat", STATUS="REPLACE",ACTION="WRITE")

             do i = 1,nx
                 WRITE(100,*) (T(i,j),j=1,end_idx-start_idx+1)
             end do
end if


 call mpi_finalize(ierr)
end program main


