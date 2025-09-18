program main
use mpi
implicit none
integer::i,x_min,n,n1,itr,j,ierr,my_rank,num_procs,start_idx,end_idx,p
real :: dt,dx,error,error_max,t,t_max,ghost_left,ghost_right
real,allocatable :: x(:),u(:),u_old(:),u_exact(:),x_position(:),u_exact_final(:),u_final(:)


call mpi_init(ierr)
call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierr)
p=num_procs

x_min=0
dx=0.002
n=int(2.0/dx)+1
print*,'n',n
start_idx = 1 + (my_rank * n) / p
end_idx = (my_rank + 1) * n / p

  allocate(x(end_idx - start_idx + 1))
  allocate(u(end_idx - start_idx + 1))
  allocate(u_old(end_idx - start_idx + 1))
  allocate(u_exact(end_idx - start_idx + 1))
  allocate(x_position(n+1))
  allocate(u_exact_final(n+1),u_final(n+1))
!print*,'harish'
do i = start_idx, end_idx
    x(i - start_idx + 1) = x_min + (i - 1) * dx
!    u(i - start_idx + 1) = x(i - start_idx + 1) * tan(x(i - start_idx + 1))
end do
do i = start_idx, end_idx
    if((x(i - start_idx + 1) > 0.0).and. (x(i - start_idx + 1).le.0.5)) then
      u(i - start_idx + 1) =sin(4*3.14*x(i - start_idx + 1) )
      else
       u(i- start_idx + 1)=0
    end if
end do
!  if (my_rank==0) then
!    write(*,*)'the start is ',start_idx
!  end if
!  if (my_rank==1) then
!    write(*,*)'the u for rank 1 is',start_idx
!  end if

dt=0.0001


itr=0
t_max=0.5
t=dt
  do i = start_idx, end_idx
     if(((x(i - start_idx + 1)-t_max)> 0.0).and. ((x(i - start_idx + 1)-t_max).le.0.5)) then
        u_exact(i - start_idx + 1) =sin(4*3.14*(x(i - start_idx + 1)-t_max) )
     else
         u_exact(i- start_idx + 1)=0.0
     end if
  end do
!print*,'ha yes?'




do while(t<=t_max)
 !print*,'ha yes?2'
!u_old=u
!if (my_rank > 0) then
!    call MPI_Send(u(1), 1, MPI_REAL, my_rank - 1, 0, MPI_COMM_WORLD, ierr)
!    call MPI_Recv(ghost_left, 1, MPI_REAL, my_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
!  end if

  if (my_rank < num_procs - 1) then
    call MPI_Send(u(end_idx - start_idx + 1), 1, MPI_REAL, my_rank+1, 0, MPI_COMM_WORLD, ierr)
   ! call MPI_Recv(ghost_right, 1, MPI_REAL, my_rank +1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  end if
  if (my_rank>0) then
    call MPI_Recv(ghost_right, 1, MPI_REAL, my_rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  end if
  if (my_rank > 0) then
    u(0) = ghost_right
  end if
!print*,'please '
!  if (my_rank < num_procs - 1) then
!    u(end_idx - start_idx + 2) = ghost_left
!  end if
!if (my_rank==0) then
!  write(*,*),'this is',u(end_idx - start_idx + 1)
! end if
if(t==0.499770641) then
if (my_rank==1) then
  write(*,*)'the u of 0 now is ',u(200)
  end if
end if
u_old=u
u_old(0)=u(0)
if (t==0.499770641) then
if(my_rank==1) then
  write(*,*)'the u old of 0 is',u_old(200)
end if
end if
!print*,'hello'
itr=itr+1
!u(1)=0
! do i=1,size(u)
    
! end do

if (my_rank==0) then
u(1)=0
 do i=2,end_idx - start_idx + 1
   u(i)=u_old(i)-(dt/dx)*(u_old(i)-u_old(i-1))
 end do
 

else if(my_rank==num_procs-1) then
do i=1,end_idx - start_idx+1 
     u(i)=u_old(i)-(dt/dx)*(u_old(i)-u_old(i-1))
 end do
 u(n)=0


else
  do i=1,end_idx - start_idx + 1
       u(i)=u_old(i)-(dt/dx)*(u_old(i)-u_old(i-1))
  end do   
end if

!do i = start_idx, end_idx
!    if(((x(i - start_idx + 1)-t)> 0.0).and. ((x(i - start_idx + 1)-t).le.0.5)) then
!      u_exact(i - start_idx + 1) =sin(4*3.14*(x(i - start_idx + 1)-t) )
!      else
!       u_exact(i- start_idx + 1)=0.0
!       end if 1000
!  end do
! u(n)=0
!print*,'the t is',t

t=t+dt
print*,'itre is',itr
end do

print*,'the n is',n/p
call mpi_gather(x,n/p,MPI_REAL,x_position,n/p,MPI_REAL,0,mpi_comm_world,ierr)
call mpi_gather(u,n/p,MPI_REAL,u_final,n/p,MPI_REAL,0,mpi_comm_world,ierr)
call mpi_gather(u_exact,n/p,MPI_REAL,u_exact_final,n/p,MPI_REAL,0,mpi_comm_world,ierr)
print*,'n/p is',n/p
!if (my_rank==0) then
!print*,'all_position',x_position
!end if
!print *, "x        Numerical Solution        Exact Solution"
!if (my_rank==1) then

 ! print *, "x        Numerical Solution     "
  ! do i=1,end_idx-start_idx+1
      !  print *, x(i) , u(i)
   ! end do !, u_exact(i)
   ! print *, "Number of iterations:", itr

! end if
! if (my_rank==1) then
!    OPEN(UNIT=100, FILE="u_assgn_2_parallel.dat", STATUS="REPLACE", ACTION="WRITE")
!
   !     do i = 1,end_idx - start_idx + 1
   !            WRITE(100,*) u(i)
  !       end do
 !   OPEN(UNIT=100, FILE="u_exact_asign_2_parallel.dat", STATUS="REPLACE", ACTION="WRITE")
!
 !        do i = 1,end_idx - start_idx + 1
 !              WRITE(100,*) u_exact(i)
 !        end do
!    OPEN(UNIT=100, FILE="x_i_parallel.dat", STATUS="REPLACE", ACTION="WRITE")
!         do i= 1,end_idx - start_idx + 1
!               WRITE(100,*) x(i)
!         end do
!   end if

   if (my_rank==0) then
    OPEN(UNIT=100, FILE="u_assgn_2_parallel_1.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1,1000
               WRITE(100,*) u_final(i)
         end do
    OPEN(UNIT=100, FILE="u_exact_asign_2_parallel_1.dat", STATUS="REPLACE", ACTION="WRITE")

         do i = 1,1000
               WRITE(100,*) u_exact_final(i)
         end do
    OPEN(UNIT=100, FILE="x_i_parallel_1.dat", STATUS="REPLACE", ACTION="WRITE")
         do i= 1,1000
               WRITE(100,*) x_position(i)
         end do
   end if
      
call mpi_finalize(ierr)
end program main

