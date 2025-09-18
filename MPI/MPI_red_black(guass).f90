program main
use mpi
implicit none
integer :: itr,n,i,j,numprocs,my_rank,start_idx,end_idx,ierr,p
 double precision ::t1,t2, error,error_max,del,status,tol,diff,x_min,diffw,y_min,error1
 double precision ,allocatable :: x(:),y(:),phi(:,:),phi_old(:,:),q(:,:),phi_c(:),phi_d(:),phi_z(:,:),y_positionf(:),phi_check(:)

    call mpi_init(ierr)
    call mpi_comm_size(mpi_comm_world,numprocs,ierr)
    call mpi_comm_rank(mpi_comm_world,my_rank,ierr)
t1 = MPI_Wtime()
p=numprocs
del=0.01
n=int(2.0/del)+1
print*,'the n is',n
y_min=-1.0
start_idx = 1 + (my_rank * n) / p
end_idx = (my_rank + 1) * n / p

allocate(x(n+1),phi_c(n-1),phi_d(n-1))
allocate(y(end_idx - start_idx + 1))
allocate(y_positionf(n))

do i = start_idx, end_idx
    y(i - start_idx + 1) = y_min + (i - 1) * del
end do

x(1)=-1.0
do i=2,n+1
    x(i)=x(i-1)+del
end do

 
allocate(phi(n+1,end_idx-start_idx+1))
allocate(phi_old(n+1, end_idx-start_idx+1))
allocate(q(n+1, end_idx-start_idx+1))
allocate(phi_z(n+1,end_idx-start_idx+2))


do i=1,n+1
  do j=start_idx,end_idx
         phi(i,(j- start_idx + 1))=0
         phi_old(i,(j- start_idx + 1))=0
         q(i,(j- start_idx + 1))=0
       end do
  end do
  do i=1,n+1
    do j=start_idx,end_idx+1
       phi_z(i,(j- start_idx + 1))=0
    end do
  end do

if(my_rank==0) then
do i=1,n+1

      phi(i,1)=sin(2*3.14*x(i))

end do
end if
if (my_rank==0) then

OPEN(UNIT=100, FILE="phi_check.dat", STATUS="REPLACE", ACTION="WRITE")
        do i=1,n+1

                 WRITE(100,*)( phi(i,j),j = 1,end_idx-start_idx+1)

        end do

end if

   do i=start_idx,end_idx
     phi(1,(i- start_idx + 1))=0
     phi(n+1,(i- start_idx + 1))=0
   end do
do  i=1,n+1

   do j=start_idx,end_idx
    
       
        q(i,(j- start_idx + 1))= (y(j-start_idx+1)**2)+(x(i)**2)
   
   end do

end do



error_max=0.0001
diff=1
itr=0

do while (diff>error_max)
    diffw=0
    itr=itr+1
    print*,itr
     if (my_rank<numprocs-1) then
     
       
call MPI_Send(phi(2:n,end_idx-start_idx+1), n-1, mpi_double_precision, my_rank+1, 0, mpi_comm_world, ierr)
call MPI_Recv(phi_d, n-1, mpi_double_precision, my_rank+1, 0, mpi_comm_world, MPI_STATUS_IGNORE, ierr)

     end if
     if(my_rank>0) then
        
       call mpi_Send(phi(2:n,1), n-1, mpi_double_precision, my_rank-1, 0, mpi_comm_world, ierr)
       call mpi_Recv(phi_c, n-1, mpi_double_precision, my_rank-1, 0, mpi_comm_world, MPI_STATUS_IGNORE, ierr)

     end if

  if(my_rank==0) then
  do i=2,n
     do j=2,end_idx- start_idx
          if(mod((i+j),2)==0) then
                phi(i,j)=(phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+phi(i+1,j)+(q(i,j)*(del**2)))/4.0
           end if
        
       end do
   end do
  do i=2,n
    do j=2,end_idx- start_idx
         if(mod((i+j),2)==1) then
            phi(i,j)=(phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+phi(i+1,j)+(q(i,j)*(del**2)))/4.0
         end if
    end do
end do
do i=2,n



phi(i,end_idx - start_idx + 1) = &
    (phi(i,end_idx - start_idx) + phi(i+1,end_idx - start_idx + 1) + &
    phi(i-1,end_idx - start_idx + 1) + phi_d(i - 1) + &
    (q(i,end_idx - start_idx + 1) * (del**2))) / 4.0

end do

do i=2,n
    do j=2,end_idx- start_idx +1
      
           diffw=diffw+abs(phi(i,j)-phi_old(i,j))
        end do
     end do

     else if(my_rank==numprocs-1) then

     do i=2,n
       do j=2,end_idx- start_idx
          if(mod((i+j),2)==0) then
                phi(i,j)=(phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+phi(i+1,j)+(q(i,j)*(del**2)))/4.0
           end if
       
       end do
   end do
  do i=2,n
    do j=2,end_idx- start_idx
         if(mod((i+j),2)==1) then
            phi(i,j)=(phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+phi(i+1,j)+(q(i,j)*(del**2)))/4.0
         end if
    end do
end do
 

     do i=2,n
        phi(i,end_idx- start_idx+1)=(phi(i-1,end_idx- start_idx+1)+ &
        phi_z(i,end_idx- start_idx+1+1)+phi(i,end_idx- start_idx+1-1)+ &
        phi(i+1,end_idx- start_idx+1)+(q(i,end_idx- start_idx+1)*(del**2)))/4.0
     end do

     do i=2,n
         phi(i,1)=(phi_c(i-1)+phi(i+1,1)+phi(i,1+1)+phi(i-1,1)+(q(i,1)*(del**2)))/4.0
      end do
     
    do i=2,n
       phi_z (i,end_idx-start_idx+1+1)=((4*phi(i,end_idx-start_idx+1))-phi(i,end_idx-start_idx+1-1))/3.0
    end do
        diffw=0
do i=2,n
    do j=1,end_idx- start_idx +1
        diffw=diffw+abs(phi(i,j)-phi_old(i,j))
    end do
end do

    else 
  
       do i=2,n
     do j=2,end_idx- start_idx
          if(mod((i+j),2)==0) then
                phi(i,j)=(phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+phi(i+1,j)+(q(i,j)*(del**2)))/4.0
           end if
      
       end do
   end do
  do i=2,n
    do j=2,end_idx- start_idx
         if(mod((i+j),2)==1) then
            phi(i,j)=(phi(i-1,j)+phi(i,j+1)+phi(i,j-1)+phi(i+1,j)+(q(i,j)*(del**2)))/4.0
         end if
    end do
end do
     do i=2,n
         phi(i,end_idx- start_idx +1)=(phi(i-1,end_idx- start_idx +1)+ &
         phi(i+1,end_idx- start_idx +1)+phi(i,end_idx- start_idx +1-1)+ &
         phi_d(i-1)+(q(i,end_idx- start_idx +1)*(del**2)))/4.0
     end do
     do i=2,n
        phi(i,1)=(phi_c(i-1)+phi(i,1+1)+phi(i-1,1)+phi(i+1,1)+(q(i,1)*(del**2)))/4.0 
     end do

    

 diffw=0
do i=2,n
    do j=1,end_idx- start_idx +1
        diffw=diffw+abs(phi(i,j)-phi_old(i,j))
    end do
end do
end if

phi_old=phi
call mpi_allreduce(diffw,diff,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

print*,'error is',diff
 print*,'iteration is',itr


end do

call mpi_gather(y,n/p,MPI_double_precision,y_positionf,n/p,MPI_double_precision,0,mpi_comm_world,ierr)

t2 = MPI_Wtime()
print*,'the time taken is',t2-t1
print*,'the iterations are ',itr
!if (my_rank==0) then
!print*,'the y is ',y_positionf
! OPEN(UNIT=100, FILE="phi_check_f.dat", STATUS="REPLACE", ACTION="WRITE")
!do i=1,200
 
 !                 WRITE(100,*) phi_check(i)
 
 !       end do
 
!end i
if (my_rank==p/2) then
!print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_jacobi.dat", STATUS="REPLACE", ACTION="WRITE")
do i=1,n+1
 
                  WRITE(100,*) phi(i,1)
 
        end do
 
end if
!end if

if (my_rank==0) then
!print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_f0.dat", STATUS="REPLACE", ACTION="WRITE")
do j = 1,end_idx-start_idx+1
 
                  WRITE(100,*) phi(101,j)
 
        end do
 
!OPEN(UNIT=100, FILE="phi_check_f.dat", STATUS="REPLACE", ACTION="WRITE")
! do j = 1,n

 !                  WRITE(100,*) phi(101,j)

 !        end do

 end if


!print*,'harish come on come on'
if (my_rank==1) then
 !print*,'the y is ',y_positionf
  OPEN(UNIT=100, FILE="phi_check_f1.dat", STATUS="REPLACE", ACTION="WRITE")
 do j = 1,end_idx-start_idx+1

                   WRITE(100,*) phi(101,j)

         end do

 end if
if (my_rank==2) then
 !print*,'the y is ',y_positionf
  OPEN(UNIT=100, FILE="phi_check_f2.dat", STATUS="REPLACE", ACTION="WRITE")
 do j = 1,end_idx-start_idx+1
  
                   WRITE(100,*) phi(101,j)
  
         end do
 
 end if
if (my_rank==3) then
 !print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_f3.dat", STATUS="REPLACE", ACTION="WRITE")

do j = 1,end_idx-start_idx+1
  
     WRITE(100,*) phi(101,j)

end do
 
end if
if (my_rank==4) then
 !print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_f4.dat", STATUS="REPLACE", ACTION="WRITE")

do j = 1,end_idx-start_idx+1

     WRITE(100,*) phi(101,j)

end do

end if
if (my_rank==5) then
 !print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_f5.dat", STATUS="REPLACE", ACTION="WRITE")

do j = 1,end_idx-start_idx+1

     WRITE(100,*) phi(101,j)

end do

end if
if (my_rank==6) then
 !print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_f6.dat", STATUS="REPLACE", ACTION="WRITE")

do j = 1,end_idx-start_idx+1

     WRITE(100,*) phi(101,j)

end do

end if
!print*,'yeah its ok'
if (my_rank==7) then
 !print*,'the y is ',y_positionf
 OPEN(UNIT=100, FILE="phi_check_f7.dat", STATUS="REPLACE", ACTION="WRITE")

do j = 1,end_idx-start_idx+1

     WRITE(100,*) phi(101,j)

end do

end if

!print*,'come on boy'
call mpi_finalize(ierr)
end program main   

