include "/opt/intel/mkl/include/mkl_dfti.f90"

subroutine fft2d
  
    use MKL_DFTI
    use variables
!$use omp_lib
    implicit none
    
    integer iwrit,div,data_set,num_data,L2,nfft,div_start,div_end,qi,iloc
    integer::Status,LEN
    TYPE(DFTI_DESCRIPTOR),POINTER :: hand
    
    real(8) dt,fsmach,SR_St,Lref,Cref,Tfft,Fac,Rez,Imz,Tref

    
    
    !real(8),allocatable :: strauhal(:),r(:),r_PSD(:),r_PSD_total(:),wsave(:)
    
    real(8),allocatable :: psd_n_dm(:,:) 
    
    character error_message*(DFTI_MAX_MESSAGE_LENGTH)
    
 
    allocate(fft_rx(nx,ny))
    allocate(psd_n_dm(nx,ny))
    
    !pi = 3.141592d0
!    pi = 4.0d0 * atan(1.0d0)

!-----------------------2DFFT start--------------------------------- 

    !----------------------substitute data to rx----------------
    
    fft_rx(:,:) = un_box(:,:) 

    !------------------------FFT for x ---------------------------
    write(*,*)"FFT for x direction start"
 
    div = nx
    allocate(rx(nx))
    rx = 0.0d0
    !---FFT preparation
    Status = DftiCreateDescriptor(hand,DFTI_DOUBLE,DFTI_COMPLEX,1,div)
    write(*,*)"DftiCreateDescriptor",Status
    Status = DftiCommitDescriptor(hand)
    write(*,*)"DftiCommitDescriptor",Status
    !Status = DftiSetValue(hand,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
    !write(*,*)"DftiSetValue",Status
    
    !--------------------




    do ii=1,ny
       
       !---FFT---------------
       
       rx(:) = fft_rx(:,ii)
       Status = DftiComputeForward (hand,rx)
  
       !write(*,*)'status',Status
       !write(error_message, '(A)')DftiErrorMessage(Status)
       !write(*,*)error_message
       !
       !------------------------
       
       !---contain FFT result 
  !     fft_rx(ii,:) = rx(:)  
       
       fft_rx(:,ii) = rx(:)  
       !---------------------
       
       
    end do

    deallocate(rx)

    
    write(*,*)"FFT end"
    
    write(*,*)"FFT for y direction start"

    div = ny

    allocate(rx(ny))
    rx = 0.0d0
    !---FFT preparation
    Status = DftiCreateDescriptor(hand,DFTI_DOUBLE,DFTI_COMPLEX,1,div)
    write(*,*)"DftiCreateDescriptor",Status
    Status = DftiCommitDescriptor(hand)
    write(*,*)"DftiCommitDescriptor",Status
    !Status = DftiSetValue(hand,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
    !write(*,*)"DftiSetValue",Status
    
    !--------------------

    
    do ii=1,nx
       
       !---FFT---------------
       rx(:) = fft_rx(ii,:)     
       Status = DftiComputeForward(hand,rx)
       

       fft_rx(ii,:) = rx(:)
       write(*,*)'status',Status
       write(error_message, '(A)')DftiErrorMessage(Status)
       write(*,*)error_message
       !


       !------------------------

       !---contain FFT result 
       
        fft_rx(ii,:) = rx(:)  
       ! fft_rx(:,ii) = rx(:)  
      
       !---------------------
       
       
    end do


    deallocate(rx)    
    !--------------------FFT end---------------------------------------------


    !-------------------make 2DFFT graph-------------------------------------
    

    !-------------------make freqency&wavenumber spectra---------------------

      call dftf(nx,kx,hx); call dftf(ny,ky,hy)

    !-------------------------------------------------------------------------


    !------------------sort FFT result in frequency order--------------------

      call dfts2d(nx,ny,kx,ky,fft_rx,mkx,mky,mz) !sort frequency

    !------------------------------------------------------------------------

    !make Power Spectra Density

    
    psd_n_dm(:,:) = cdabs( fft_rx(:,:) ) * cdabs( fft_rx(:,:) ) 

    psd_n_dm(:,:) = psd_n_dm(:,:) / ( len_x * len_y  )
    

    
    !------------------------------------------------------------------------
    
   
    !---------------output the file with pm3d format (for gnuplot)------------

    open(21,file="FFT2d_pm3d.txt")
    do i=0,Nx-1
       do j=0,Ny-1
          write(21,'(4e20.10e3)')mkx(i),mky(j),cdabs(mz(i,j))
          !write(21,*)i,j,cdabs(mz(i,j))
          
       enddo
       write(21,*)
    enddo
    close(21)
    
    !------------------------------------------------------------------------



end subroutine






subroutine dftf(N,f,h)
  implicit none
  integer,intent(in)::N
  double precision,intent(in)::h
  double precision,intent(out)::f(0:N-1)

  integer::i
  double precision::mf(0:N-1)
 
  if(mod(N,2).eq.0)then
     do i=0,N-1
        mf(i)=(dble(i+1)-dble(N)*0.5d0)/(dble(N)*h)
     enddo

     do i=0,N-1
        if(i.le.N/2-2)then
           f(i+N/2+1)=mf(i)
        else
           f(i-N/2+1)=mf(i)
        endif
     enddo
  else
     do i=0,N-1
        mf(i)=(dble(i)-dble(N-1)*0.5d0)/(dble(N)*h)
     enddo

     do i=0,N-1
        if(i.le.(N-1)/2-1)then
           f(i+(N-1)/2+1)=mf(i)
        else
           f(i-(N-1)/2)=mf(i)
        endif
     enddo
  endif
 
  return
end subroutine dftf



subroutine dfts2d(Nx,Ny,fx,fy,z,mfx,mfy,mz)
  implicit none
  integer,intent(in)::Nx,Ny
  double precision,intent(in)::fx(0:Nx-1),fy(0:Ny-1)
  complex(kind(0d0))::z(0:Nx-1,0:Ny-1)
  double precision,intent(out)::mfx(0:Nx-1),mfy(0:Ny-1)
  complex(kind(0d0)),intent(out)::mz(0:Nx-1,0:Ny-1)
  complex(kind(0d0))::mmz(0:Nx-1,0:Ny-1)
  integer::i,j,k,l  


  if(mod(Ny,2).eq.0)then
     do i=0,Ny-1
        if(i.le.Ny/2)then
           j=i+Ny/2-1
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        else
           j=i-Ny/2-1
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        endif
     enddo
  else
     do i=0,Ny-1
        if(i.le.(Ny-1)/2)then
           j=i+(Ny-1)/2
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        else
           j=i-(Ny-1)/2-1
           mfy(j)=fy(i)
           mz(0:Nx-1,j)=z(0:Nx-1,i)
        endif
     enddo
  endif

  mmz=mz

  if(mod(Nx,2).eq.0)then
     do k=0,Nx-1
        if(k.le.Nx/2)then
           l=k+Nx/2-1
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        else
           l=k-Nx/2-1
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        endif
     enddo
  else
     do k=0,Nx-1
        if(k.le.(Nx-1)/2)then
           l=k+(Nx-1)/2
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        else
           l=k-(Nx-1)/2-1
           mfx(l)=fx(k)
           mz(l,0:Ny-1)=mmz(k,0:Ny-1)
        endif
     enddo
  endif  

  return
end subroutine dfts2d
