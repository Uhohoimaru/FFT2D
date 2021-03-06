

 module m_bmp
      implicit none
      type :: t_bmp_file_header
        sequence  
        integer(2) :: bfType = transfer('BM', 0_2, 1) ! BitMap
        integer(4) :: bfSize          ! file size in bytes
        integer(2) :: bfReserved1 = 0 ! always 0
        integer(2) :: bfReserved2 = 0 ! always 0
        integer(4) :: bfOffBits
      end type t_bmp_file_header
      !
      type :: t_bmp_info_header
        sequence
        integer(4) :: biSize     = Z'28' ! size of bmp_info_header ; 40bytes 
        integer(4) :: biWidth
        integer(4) :: biHeight
        integer(2) :: biPlanes   = 1 ! always 1
        integer(2) :: biBitCount
        integer(4) :: biCompression = 0 !0:nocompression,1:8bitRLE,2:4bitRLE,3:bitfield
        integer(4) :: biSizeImage
        integer(4) :: biXPelsPerMeter = 2 ! 96dpi
        integer(4) :: biYPelsPerMeter = 2 ! 96dpi 
        integer(4) :: biClrUsed      = 0
        integer(4) :: biClrImportant = 0 
      end type t_bmp_info_header
      !
      type :: t_rgb
        sequence
        character :: b, g, r
      end type t_rgb 
      type :: t_bmp
        type (t_rgb), allocatable :: rgb(:, :)
      contains 
        procedure :: init => init_bmp
        procedure :: wr => wr_bmp
        procedure :: rd => rd_bmp
      end type  
    contains   
      subroutine init_bmp(bmp, nx, ny)
        class (t_bmp), intent(in out) :: bmp
        integer, intent(in) :: nx, ny
        allocate(bmp%rgb(nx, ny))
      end subroutine init_bmp

      subroutine wr_bmp(bmp, fn)
        class (t_bmp), intent(in) :: bmp
        character (len = *), intent(in) :: fn
        type (t_bmp_file_header) :: bmp_file_header
        type (t_bmp_info_header) :: bmp_info_header
        integer :: i, j, k, l
        associate(nx => size(bmp%rgb, 1), ny => size(bmp%rgb, 2))
          bmp_file_header%bfSize      = 14 + 40 + 0 + (3 * nx + mod(nx, 4)) * ny
          bmp_file_header%bfOffBits   = 14 + 40
          bmp_info_header%biWidth     = nx
          bmp_info_header%biHeight    = ny
          bmp_info_header%biBitCount  = 24 
          bmp_info_header%biSizeImage = (3 * nx + mod(nx, 4)) * ny
          open(9, file = fn//'.bmp', access = 'stream', status = 'unknown')
          write(9) bmp_file_header
          write(9) bmp_info_header
          write(9) (bmp%rgb(:, i), (achar(0), j = 1, mod(nx, 4)), i = 1, ny)
          close(9)
        end associate
      end subroutine wr_bmp

      subroutine rd_bmp(bmp, fn)
        class (t_bmp), intent(out) :: bmp
        character (len = *), intent(in) :: fn
        type (t_bmp_file_header) :: bmp_file_header
        type (t_bmp_info_header) :: bmp_info_header
        integer :: i, j, k,l
        character :: dummy,num_r, num_g, num_b
        real(8),allocatable :: num_gray(:,:)
 

        associate(nx => bmp_info_header%biWidth, ny => bmp_info_header%biHeight)
          open(10, file = fn//'.bmp', access = 'stream', status = 'old')
          read(10) bmp_file_header
          read(10) bmp_info_header
          allocate(bmp%rgb(nx, ny))
          read(10) (bmp%rgb(:, i), (dummy, j = 1, mod(nx, 4)), i = 1, ny)
           close(10)

        allocate(num_gray(nx,ny))
        
           write(*,*)"a",iachar("a")
           

          open(150,file="Lenna_test.txt")
          
          do j=1,nx
             do l=1,ny
                write(150,*)iachar(bmp%rgb(j,l)%r),iachar(bmp%rgb(j,l)%g),iachar(bmp%rgb(j,l)%b)

                !!!!Gray scale!!!(ImageMagick)

                num_gray(j,l) = ( dble( iachar(bmp%rgb(j,l)%r) )*0.298839 + dble( iachar(bmp%rgb(j,l)%g) )*0.586811 + dble( iachar(bmp%rgb(j,l)%b) )*0.114350 ) /256
                
                !!!!!!!!!!!!!!!!!
                
             end do
          end do
          
          open(20,file="./input_size.txt")
          write(20,*)nx,ny
          close(20)

          open(200,file="./Lenna_grayscale_pm3d.txt")
          
          do l=1,ny
             do j=1,nx
   
                   write(200,*)j,l,num_gray(j,l)
                
                end do
                write(200,*)

             end do

             close(200)
            
        end associate  


      end subroutine rd_bmp
    end module m_bmp







    subroutine readdata_forbmp
      use m_bmp
      use variables
      implicit none
      type (t_bmp) :: pic0, pic1
      call pic0%rd('Lenna')
      allocate(pic1%rgb, mold = pic0%rgb)
      pic1%rgb = pic0%rgb(size(pic0%rgb, 1):1:-1, :) ! reverse
      call pic1%wr('reverse')

    end subroutine  readdata_forbmp





