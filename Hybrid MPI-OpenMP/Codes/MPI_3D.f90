program main
use mpi
implicit none

    ! Input parameters
    integer ::itr, nx, ny,nz,nt,i,j,k,my_rank,num_procs,p,start_idx,end_idx,ierr
    real ::L, dt, alpha,time,dx,dy,dz,k1,k2,y_min,k3
    double precision:: t1,t2
    real, allocatable :: x(:),y(:),z(:),T_c(:,:),T_d(:,:)
    real, dimension(:,:,:), allocatable :: T, Told,all_T

  
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
    L=100
    nx = 64
    ny = 64
    nz = 64
    dt = 0.02
    alpha = 1.4
    nt = 2000
    dx=L/nx
    dy=L/ny
    dz=L/nz
start_idx = 1 + (my_rank * ny) / p
end_idx = (my_rank + 1) * ny / p

allocate(y(end_idx - start_idx + 1),x(nx),z(nz),T_c(nx,nz),T_d(nx,nz))
   x(1)=0
   y_min=0
   z(1)=0
!    y(1)=0
     do i=2,nx

       x(i)=x(i-1)+dx

     end do

     do i = start_idx, end_idx
          y(i - start_idx + 1) = y_min + (i - 1) * dy
     end do

     do i= 2,nz
        
        z(i)=z(i-1)+dz
 
     end do

!     do i=2,ny

 !      y(i)=y(i-1)+dy

  !   end do
     if(my_rank==0) then
        print*,'harish',start_idx
        print*,'harish',end_idx

     end if
    allocate(T(nx,end_idx-start_idx+1,nz))
    allocate(Told(nx,end_idx-start_idx+1,nz))
    allocate(all_T(nx, ny,nz))
   print*,'harish22'

   T = 300.0
    
    do i = start_idx,end_idx
        T(1,i-start_idx+1,:) = 600.0
        T(ny,i-start_idx+1,:) = 900.0
       ! T(:,i-start_idx+1,1)=700
       ! T(:,i-start_idx+1,nz)=500

    end do

if (my_rank==0) then
    do i = 1, nx
        T(i,1,:) = 400.0
      ! T(j,nx) = 800.0
    end do
end if

if (my_rank==num_procs-1) then
  do i=1,nx
      T(i,end_idx-start_idx+1,:)=800.0
  end do
end if
    
    do i = start_idx,end_idx
       T(:,i-start_idx+1,1)=700
       T(:,i-start_idx+1,nz)=500
    end do

if (my_rank==1) then
        do i = 1,nx
            WRITE(*,*) 'the temp is',(T(i,j,1),j=1,end_idx-start_idx+1)
         end do
    end if
    Told=T
     if (my_rank==2) then
        do i = 1,nx
           ! WRITE(*,*) 'the temparature old in rannk 0 is',(Told(i,j),j=1,end_idx-start_idx+1)
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
    k3 = alpha * (dt / (dz**2))
   ! CFL = k1 + k2
   time=0;
   itr=0;
    ! Time loop
do while (time< nt)
itr=itr+1
 if (my_rank==1) then

!print*,'the itr is',itr
end if
     if (my_rank < num_procs-1) then
        call MPI_Send(T(1:nx,end_idx-start_idx+1,:), nx*nz, mpi_real, my_rank+1, 0, mpi_comm_world, ierr)
        call MPI_Recv(T_d, nx*nz, mpi_real, my_rank+1, 0, mpi_comm_world, MPI_STATUS_IGNORE, ierr)

     end if
     if(my_rank>0) then

           call mpi_Send(T(1:nx,1,:), nx*nz, mpi_real, my_rank-1, 0, mpi_comm_world, ierr)
           call mpi_Recv(T_c, nx*nz, mpi_real, my_rank-1, 0, mpi_comm_world, MPI_STATUS_IGNORE, ierr)

     end if

       if(my_rank==0) then

          do i = 2, nx-1
            do j = start_idx+1, end_idx-1
	       do k=2,nz-1
                T(i,j-start_idx+1,k) = Told(i,j-start_idx+1,k) + (k1 * (Told(i+1,j-start_idx+1,k) - 2.0*Told(i,j-start_idx+1,k)&
                                    + Told(i-1,j-start_idx+1,k))) &
                                  + (k2 * (Told(i,j-start_idx+1-1,k) - 2.0*Told(i,j-start_idx+1,k) + Told(i,j-start_idx+1+1,k)))&
				  + k3 * (Told(i,j-start_idx+1,k+1) - 2.0*Told(i,j-start_idx+1,k) + Told(i,j-start_idx+1,k-1))
            end do
          end do
        end do
         do i = 2, nx-1
            do k=2,nz-1

                T(i,end_idx-start_idx+1,k) = Told(i,end_idx-start_idx+1,k) + (k1 * (Told(i+1,end_idx-start_idx+1,k)&
                                   - 2.0*Told(i,end_idx-start_idx+1,k) + Told(i-1,end_idx-start_idx+1,k))) &
                                  + (k2 * (Told(i,end_idx-start_idx+1-1,k) - 2.0*Told(i,end_idx-start_idx+1,k) + T_d(i,k)))&
				  + k3 * (Told(i,end_idx-start_idx+1,k+1) - 2.0*Told(i,end_idx-start_idx+1,k) + Told(i,end_idx-start_idx+1,k-1))
            end do
	    end do
      else if (my_rank==num_procs-1) then
          do i = 2, nx-1
            do j = start_idx+1,end_idx-1
	       do k=2,nz-1
                T(i,j-start_idx+1,k) = Told(i,j-start_idx+1,k) + (k1 * (Told(i+1,j-start_idx+1,k) - 2.0*Told(i,j-start_idx+1,k)&
                                  + Told(i-1,j-start_idx+1,k))) &
                                  + (k2 * (Told(i,j-start_idx+1-1,k) - 2.0*Told(i,j-start_idx+1,k) + Told(i,j-start_idx+1+1,k)))&
				  + k3 * (Told(i,j-start_idx+1,k+1) - 2.0*Told(i,j-start_idx+1,k) + Told(i,j-start_idx+1,k-1))
            end do
        end do
         end do
         do i=2,nx-1
	   do k=2,nz-1
                
               T(i,1,k)=Told(i,1,k) + (k1 * (Told(i+1,1,k) - 2.0*Told(i,1,k) + Told(i-1,1,k))) &
                                  + (k2 * (T_c(i,k) - 2.0*Told(i,1,k) + Told(i,1+1,k)))&
				 + k3 * (Told(i,1,k+1) - 2.0*Told(i,1,k) + Told(i,1,k-1))
         end do
	 end do
      else

         do i = 2, nx-1
            do j=start_idx+1,end_idx-1
	       do k=2,nz-1
                T(i,j-start_idx+1,k) =  Told(i,j-start_idx+1,k) + (k1 * (Told(i+1,j-start_idx+1,k)&
                                    - 2.0*Told(i,j-start_idx+1,k) + Told(i-1,j-start_idx+1,k))) &
                                    + (k2 * (Told(i,j-start_idx+1-1,k)- 2.0*Told(i,j-start_idx+1,k) + Told(i,j-start_idx+1+1,k)))&
				    + k3 * (Told(i,j-start_idx+1,k+1) - 2.0*Told(i,j-start_idx+1,k) + Told(i,j-start_idx+1,k-1))
            end do
        end do
        end do
         do i=2,nx-1
            do k=2,nz-1
             T(i,1,k)=Told(i,1,k) + (k1 * (Told(i+1,1,k) - 2.0*Told(i,1,k) + Told(i-1,1,k))) &
                             + (k2 * (T_c(i,k) - 2.0*Told(i,1,k) + Told(i,1+1,k)))&
			     + k3 * (Told(i,1,k+1) - 2.0*Told(i,1,k) + Told(i,1,k-1))
         end do
         end do
        do i=2,nx-1
            do k=2,nz-1
              T(i,end_idx-start_idx+1,k) = Told(i,end_idx-start_idx+1,k) + (k1 * (Told(i+1,end_idx-start_idx+1,k)&
                                       - 2.0*Told(i,end_idx-start_idx+1,k) + Told(i-1,end_idx-start_idx+1,k))) &
                                       + (k2 * (Told(i,end_idx-start_idx+1-1,k) - 2.0*Told(i,end_idx-start_idx+1,k) + T_d(i,k)))&
                                       + k3 * (Told(i,end_idx-start_idx+1,k+1) - 2.0*Told(i,end_idx-start_idx+1,k)&
				       + Told(i,end_idx-start_idx+1,k-1))
         end do
         end do
    end if

        Told = T
        time=time+dt

end do
t2=mpi_wtime()
print*,'the time taken is',t2-t1
        do i = 1,nx
  !          WRITE(*,*) 'the temparature final in rank 2 is',(T(i,j),j=1,end_idx-start_idx+1)
         end do
 if (my_rank==0) then
       do i=1,nx

!            write(*,*) (T(i,j,2),j=1,end_idx-start_idx+1)
   
         end do
	 end if
     call MPI_Allgather(T, nx*ny*nz/p, mpi_real, all_T, nx*ny*nz/p, mpi_real, MPI_COMM_WORLD, ierr)

 print*, 'the time is',time
OPEN(UNIT=100, FILE="time_mpi_3d.dat", STATUS="REPLACE", ACTION="WRITE")
WRITE(100,*) t2-t1

OPEN(UNIT=100, FILE="num_procs_mpi_3d.dat", STATUS="REPLACE", ACTION="WRITE")
WRITE(100,*) p

if (my_rank==(p/2)-1) then
OPEN(UNIT=100, FILE="temperature_mpi_3d.dat", STATUS="REPLACE", ACTION="WRITE")

 do i=1,nx

   write(100,*) (T(i,end_idx-start_idx+1,k),k=1,nz)

end do
end if
print*,'itr is',itr
!OPEN(UNIT=100, FILE="x.dat", STATUS="REPLACE", ACTION="WRITE")

!         do i = 1,nx
 !              WRITE(100,*) x(i)
  !       end do
!OPEN(UNIT=100, FILE="y.dat", STATUS="REPLACE", ACTION="WRITE")
  
   !        do i = 1,nx
    !             WRITE(100,*) y(i)
     !      end do


!OPEN(UNIT=100, FILE="temperature_mpi_3d.dat", STATUS="REPLACE", ACTION="WRITE")
 
    !         do i = 1,nx
     !            WRITE(100,*) (T(i,2,k),k=1,nz)
      !       end do
 call mpi_finalize(ierr)
end program main


