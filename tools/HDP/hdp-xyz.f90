!clustering by fast search-and-find of density peak
!calculate density from atomic position
!written by Y. S
module global
   integer :: ndata
   real :: dcut,dcut_2
   real,parameter :: zero_point_cut=0.5
   real,allocatable :: delta(:),x(:,:),rou(:),pool(:),distance2(:)
end module global

program main
   use global
   implicit none
   integer i,j,k,temp
   real rou_max,delta_max
   real,allocatable :: index_sequence(:)
   integer,allocatable :: list(:)
   read(*,*) dcut
   dcut_2=dcut**2
   call inp
   call rou_cal
   call delta_cal
   allocate(list(ndata),index_sequence(ndata))
   open(99,file='density-delta.dat',status='unknown')
   rou_max=maxval(rou)
   delta_max=maxval(delta)
   do i=1,ndata
      index_sequence(i)=(rou(i)/rou_max)**2+(delta(i)/delta_max)**2
   enddo
   do i=1,ndata
      list(i)=i
   enddo
   do i=1,ndata-1
      do j=1,ndata-1-i
        ! if(index_sequence(list(j)) .lt. index_sequence(list(j+1)) )then
         if(delta(list(j)) .lt. delta(list(j+1)) )then
            temp=list(j)
            list(j)=list(j+1)
            list(j+1)=temp
         endif
      enddo
   enddo
   do k=1,ndata
      i=list(k)
      write(99,'(f10.1,f16.3,5f16.6)') rou(i),delta(i),x(i,1),x(i,2),x(i,3),sqrt(x(i,1)**2+x(i,2)**2+x(i,3)**2),index_sequence(i)
   !   write(99,*) rou(i),delta(i),x(i,1),x(i,2),x(i,3),index_sequence(i)
   enddo
   deallocate(rou,x,delta,index_sequence,list)
   close(99)
end program main

subroutine dis_cal(atom,poolsize)
   use global
   implicit none
   integer atom,poolsize,i,j,k
   real*8 rsq
   distance2=0
   do i=1,poolsize
      rsq=0
      do k=1,3
         rsq=rsq+(x(pool(i),k)-x(atom,k))**2
      enddo
      distance2(i)=rsq
   enddo
end subroutine

subroutine delta_cal
   use global
   implicit none
   integer i,j,k,max_no,nb,counter,rember(ndata),ct
   max_no=sum(maxloc(rou))
   nb=ndata-1
   allocate(pool(nb),distance2(nb))
   ct=1
   do i=1,ndata
      if(i.eq.max_no) cycle
      pool(ct)=i
      ct=ct+1
   enddo
   call dis_cal(max_no,nb)
   delta(max_no)=sqrt(maxval(distance2))
   deallocate(pool,distance2)
   print*,  'Max: ',rou(max_no),delta(max_no)

   do i=1,ndata
      if(i.eq.max_no) cycle
      counter=0
      rember=0
      do j=1,ndata
         if(rou(j).gt.rou(i))then
            counter=counter+1
            rember(counter)=j
         endif
      enddo
      nb=counter
      allocate(pool(nb),distance2(nb))
      do j=1,nb
         pool(j)=rember(j)
      enddo
      call dis_cal(i,nb)
      delta(i)=sqrt(minval(distance2))
      deallocate(pool,distance2)
   enddo
end subroutine delta_cal



subroutine rou_cal
   use global
   implicit none
   integer i,j,k
   real*8 ss
   rou=0
   do i=1,ndata
      do j=1,ndata
         ss=(x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+(x(i,3)-x(j,3))**2
         if(ss.lt.dcut_2)then
            rou(i)=rou(i)+1
         endif
      enddo
   enddo
end subroutine rou_cal

subroutine inp
   use global
   implicit none
   integer i,j,k,m1,m2,m3,ct,ngrd,nli,ntotal
   real*8 box,delta_r,d(5)
  ! real*8,allocatable :: tdensity(:,:,:)
   character*2 aname
!   open(5, access='sequential')
   open(5,file='hdp.xyz',status='old')
   read(5,*) ndata
   read(5,*)
   allocate(rou(ndata),x(ndata,3),delta(ndata))
   do i=1,ndata
      read(5,*) aname,x(i,1),x(i,2),x(i,3)
   enddo
   close(5)
end subroutine
