!--------------------------------------------
! created by: MKL sample code (cite:https://software.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-fortran/top/appendix-e-code-examples/fourier-transform-functions-code-examples/fft-code-examples.html)
! 
! rearreanged by Shota Morita, 05/10/2021
!
!------------------------------------------------

include "/opt/intel/mkl/include/mkl_dfti.f90"


subroutine fft2d
  
    use MKL_DFTI
    use variables
!$use omp_lib
    implicit none


    
    integer iwrit,div,data_set,num_data,L2,nfft,div_start,div_end,qi,iloc
    integer::Status,LEN
    TYPE(DFTI_DESCRIPTOR),POINTER :: hand
    
    real(8) dt,fsmach,un_box_ave,SR_St,Lref,Cref,Tfft,Fac,Rez,Imz,Tref,len_x,len_y

    
    
    !real(8),allocatable :: strauhal(:),r(:),r_PSD(:),r_PSD_total(:),wsave(:)
    complex(kind(0d0)),allocatable :: rx(:), rx2d(:,:) 
    !Equivalence (rx2d,rx)

    real(8),allocatable :: psd_n_dm(:,:)  
     
    character error_message*(DFTI_MAX_MESSAGE_LENGTH)
    

    INTEGER :: STRIDE(2)

    type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
    type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim2
 
    len_x = 0.01d0
    len_y = 0.01d0
    allocate(rx2d(nx,ny))
    allocate(rx(nx*ny))
    allocate(psd_n_dm(nx,ny))
     
    !pi = 3.141592d0
!    pi = 4.0d0 * atan(1.0d0)
 

    !----------------------substitute data to rx----------------
    
    rx2d(:,:) = un_box(:,:) 

         do l=1,ny
      do j=1,nx
  
      k = (l-1)*j + j

       rx(k) = rx2d(j,l) 
      
    end do
 end do
 

    !----

!-------FFT preparation start


Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,&
                                 DFTI_COMPLEX, 1, nx )
Status = DftiCreateDescriptor(Desc_Handle_Dim2, DFTI_DOUBLE,&
                                 DFTI_COMPLEX, 1, ny )

!-------FFT start------------

! perform nx one-dimensional transforms along 1st dimension
Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, ny )
Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, nx )
Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, nx )
Status = DftiCommitDescriptor( Desc_Handle_Dim1 )
Status = DftiComputeForward( Desc_Handle_Dim1, rx )


! perform ny one-dimensional transforms along 2nd dimension
Stride(1) = 0; Stride(2) = nx
Status = DftiSetValue( Desc_Handle_Dim2, DFTI_NUMBER_OF_TRANSFORMS, nx )
Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_DISTANCE, 1 )
Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_DISTANCE, 1 )
Status = DftiSetValue( Desc_Handle_Dim2, DFTI_INPUT_STRIDES, Stride )
Status = DftiSetValue( Desc_Handle_Dim2, DFTI_OUTPUT_STRIDES, Stride )
Status = DftiCommitDescriptor( Desc_Handle_Dim2 )
Status = DftiComputeForward( Desc_Handle_Dim2, rx )
Status = DftiFreeDescriptor( Desc_Handle_Dim1 )
Status = DftiFreeDescriptor( Desc_Handle_Dim2 )


!------FFT end



!00000000

   do l=1,ny
      do j=1,nx

      k = (l-1)*j + j

      rx2d(j,l) = rx(k)
      
   end do
end do



    !-------------------make 2DFFT graph-------------------------------------
    
    !make Power Spectra Density
    psd_n_dm(:,:) = cdabs( rx2d(:,:) ) * cdabs( rx2d(:,:) ) 

    psd_n_dm(:,:) = psd_n_dm(:,:) / ( len_x * len_y  )
    

    
    !------------------------------------------------------------------------

    

    open(30,file='FFT2d_pm3d.txt')
    
    do l=1,ny
       do j=1,nx
       write(30,*)dble(j)*len_x,dble(l)*len_y, psd_n_dm(j,l)
       end do
       write(30,*)
    end do
    

    
  
















end subroutine
