program main
  
  use variables
  
  implicit none
  real(8) :: un_box_ave
  

  call readdata_forbmp

!!!!read data
  
  open(20,file="./input_size.txt")
  read(20,*)nx,ny
  close(20)
  allocate(un_box(nx,ny))
  
  open(200,file="./Lenna_grayscale_pm3d.txt")
  
  do l=1,ny
     do j=1,nx
        
        read(200,*)ii,ii,un_box(j,l)
        
     end do
     read(200,*)
     
  end do
!!!!!

  write(*,*)maxval(un_box)

  write(*,*)minval(un_box)    
    !--------test
!  Nx=100; Ny=100

 
!  allocate(un_box(nx,ny))
  allocate(x(0:Nx-1),kx(0:Nx-1),mkx(0:Nx-1))
  x=0d0; kx=0d0; mkx=0d0
  allocate(y(0:Ny-1),ky(0:Ny-1),mky(0:Ny-1))
  y=0d0; ky=0d0; mky=0d0
  allocate(z(0:Nx-1,0:Ny-1),mz(0:Nx-1,0:Ny-1))
  z=dcmplx(0d0,0d0); mz=dcmplx(0d0,0d0)
  
  xmin=0
  xmax=1
  hx=(xmax-xmin)/dble(Nx)
  ymin=0
  ymax=1
  hy=(ymax-ymin)/dble(Ny)
 
  len_x = hx

  len_y = hy

  do i=0,Nx-1
     x(i)=xmin+dble(i)*hx
  enddo
  do j=0,Ny-1
     y(j)=ymin+dble(j)*hy
  enddo
!  do i=0,Nx-1
!     do j=0,Ny-1
!        un_box(i+1,j+1)=func(x(i),y(j))  
!     enddo
!  enddo

    !--------


  
  !---------substract average of sample
             un_box_ave = 0.0d0
  do l=1,nx
     do j=1,ny
     
        k=(l-1)*j + j
        un_box_ave = un_box_ave + un_box(l,j) 

     end do
  end do
  
  

  un_box_ave = un_box_ave / dble(nx*ny)
  un_box = un_box - un_box_ave
  write(*,*) un_box_ave 

      open(40,file='data_subaverage.txt')
    
    do i=1,nx
       do j=1,nx
       write(40,*)x(i-1),y(j-1), un_box(i,j)
       end do
       write(40,*)
    end do
    close(40)
  
  !-----------------------------------------
  
  call fft2d
  
  call system("sh makeplot.sh")











end program


function func(x,y)
  implicit none
  double precision::x,y
  double precision::func
 
  func=dsin(4d0*x) + dsin(3d0*y) 
 
  return
end function func




