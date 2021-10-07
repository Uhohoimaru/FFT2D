
include "/opt/intel/mkl/include/mkl_dfti.f90"
subroutine sub_fftoverlap
use MKL_DFTI
use variables
!$use omp_lib
implicit none

integer iwrit,div,data_set,num_data,L2,nfft,div_start,div_end,qi,iloc
integer::Status,LEN
TYPE(DFTI_DESCRIPTOR),POINTER :: hand

real(8) dt,fsmach,un_box_ave,SR_St,Lref,Cref,Tfft,Fac,Rez,Imz,Tref

!real(8) xre_in,xim_in,yre_out,yim_out

real(8),allocatable :: strauhal(:),r(:),r_PSD(:),r_PSD_total(:),wsave(:),xre_in(:),xim_in(:),yre_out(:),yim_out(:)

complex(kind(0d0)),allocatable :: rx(:) 


character error_message*(DFTI_MAX_MESSAGE_LENGTH)




	
	!pi = 3.141592d0
	 pi = 4.0d0 * atan(1.0d0)


	dt = deltat_nondim
	iwrit = 1
	fsmach = 0.9d0
	un_box_ave = 0.0d0
	
	SR_St = 1/(dt*dble(iwrit))

	!data separation
	div = 1024
	data_set = 18
	num_data = div*data_set

	allocate(strauhal(div))
	allocate(r(div))
	allocate(r_PSD(div))
	allocate(r_PSD_total(div))
	allocate(wsave(nbs))

	allocate(xre_in(div))
	allocate(xim_in(div))
	allocate(yre_out(div))
	allocate(yim_out(div))
	allocate(rx(div))


	L2 = div/2
	nbs = nbinst

	nfft = div
		
	Lref = 7.7163667387935219d-5
	Cref = 3.1973d2
	!Cref = uj
		
	Tref = Lref/Cref
	Tfft = dt*Tref*nfft

	Fac = 1.0d0/Tfft*Lref/Cref
	r_PSD_total = 0.0d0 

	Status = DftiCreateDescriptor(hand,DFTI_DOUBLE,DFTI_COMPLEX,1,div)
	write(*,*)"DftiCreateDescriptor",Status
  	Status = DftiCommitDescriptor(hand)
	write(*,*)"DftiCommitDescriptor",Status
	!Status = DftiSetValue(hand,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
	!write(*,*)"DftiSetValue",Status
 

	!----make strauhal number------
	do i=1,div

		!strouhal(i) = SR_St * dble(i) / dble(div)
		strauhal(i) = SR_St * dble(i-1) / dble(div)
	end do
	!----------------------------------


	do i=1,nbs

		un_box_ave = un_box_ave + un_box(i) 

	end do
	

	!---------substract average of sample
	un_box_ave = un_box_ave / dble(nbs)
	un_box = un_box - un_box_ave
	!write(*,*) un_box_ave 

	!-----------------------------------------

!------------------------FFT---------------------------

write(*,*)"FFT start"

	do qi=1,data_set

	!xre_in = 0.0d0
	!!xim_in=0.0d0
	!!yre_out=0.0d0
	!yre_out=0.0d0
	do i=1,div
	rx(i) = (0.0d0,0.0d0)
	end do
	!write(*,*)"itr =",qi,"/",data_set

		div_start = (qi-1)*div + 1
		div_end = qi*div

	!write(*,*)div_start,"~",div_end	

		!-------contain sample to r------------

		do i=1,div
		iloc = div_start + i
			!xre_in(i) = un_box(iloc)
			!xim_in(i) = 0.0d0
			rx(i) = un_box(iloc)

		end do

		!------hanning window ---------

		do i=1,div

			!xre_in(i) = xre_in(i)*(0.5d0 - 0.5d0 * cos(2.0d0*pi*dble(i)/dble(nbs)) )
			rx(i) = rx(i)*(0.5d0 - 0.5d0 * cos(2.0d0*pi*dble(i-1)/dble(div-1)) )

		end do

		!------------------------------------
	
		!---FFT---------------
		Status = DftiComputeForward (hand,rx)
		!write(*,*)'status',Status
		!write(error_message, '(A)')DftiErrorMessage(Status)
		!write(*,*)error_message
!
		!------------------------

		!----make PSD-----------
	
		do i=1,div

			rx(i) = rx(i)/ dble(div)
			!yre_out(i) = yre_out(i) / dble(div)
			!yim_out(i) = yim_out(i) / dble(div)

		end do

		!----make absolute value of r----
		
		do i=1,div

			!Rez = yre_out(i)
			!Imz = yim_out(i)
			
			!r_PSD(i) = sqrt(Rez*Rez + Imz*Imz)
			r_PSD(i) = abs(rx(i))*abs(rx(i))
	
		end do

		!-------------------------------------

		r_PSD = r_PSD / dble(SR_St/div)

		!------------------------------------

		!--summation of results
		
		r_PSD_total = r_PSD + r_PSD_total

		!-------------------------
		
		!write(*,*) maxval(r_PSD_total)
		!write(*,*) minval(r_PSD_total)
	end do

write(*,*)"FFT end"

!--------------------FFT end---------------------------------------------

!---------------------FFT overlap----------------------------------------
write(*,*)"FFT overlap start"



	do qi=1,data_set-1

	!xre_in = 0.0d0
	!xim_in=0.0d0
	!yre_out=0.0d0
	!yre_out=0.0d0
	do i=1,div
	rx(i) = (0.0d0,0.0d0)
	end do
	!write(*,*)"itr =",qi,"/",data_set-1
	
		div_start = (qi-1)*div + div/2  
		div_end = qi*div + div/2

	write(*,*)div_start,"~",div_end

		!-------contain sample to r------------

		do i=1,div
		iloc = div_start + i

			!xre_in(i) = un_box(iloc)
			!xre_in(i) = 0.0d0
			rx(i) = un_box(iloc)

		end do

		!------hanning window ---------

		do i=1,div

			!xre_in(i) = xre_in(i)*(0.5d0 - 0.5d0 * cos(2.0d0*pi*dble(i)/dble(nbs)) )
			rx(i) = rx(i)*(0.5d0 - 0.5d0 * cos(2.0d0*pi*dble(i-1)/dble(div-1)) )

		end do

		!------------------------------------
	
		!---FFT---------------
		Status = DftiComputeForward (hand,rx)
		!write(error_message, '(A)')DftiErrorMessage(Status)
		!write(*,*)error_message
		!------------------------



		!----make PSD-----------
	
		do i=1,div

			!yre_out(i) = yre_out(i) / dble(div)
			!yim_out(i) = yim_out(i) / dble(div)
			rx(i) = rx(i) / dble(div)

		end do

		!----make absolute value of r----
		
		do i=1,div

			!Rez = yre_out(i)
			!Imz = yim_out(i)
			
			!r_PSD(i) = sqrt(Rez*Rez + Imz*Imz)

			r_PSD(i) = abs(rx(i))*abs(rx(i))
	
		end do

		!-------------------------------------

			r_PSD = r_PSD / dble(SR_St/div)

	!------------------------------------

		!--summation of results
		
		r_PSD_total = r_PSD + r_PSD_total

		!-------------------------

	end do

	write(*,*)"FFT overlap end"
!------------------------FFT overlap end---------------------------

	!--------ensmable averaging-------

		r_PSD_total = r_PSD_total / dble(data_set + data_set -1)

	!---------------------------------------

	write (filename,'("./FFTresult/p/data/FFT_z=",f4.1,"r=",f4.1,".dat")')x(jpa,1,lpa),z(jpa,1,lpa)

 	open(30,file=filename,status='unknown')

	do i=2,L2
	write(30,*)Strauhal(i),r_PSD_total(i)
	end do
	
 maxPSD = 0.0d0
        do i=2,L2

           if(maxPSD < r_PSD_total(i)) then

              maxPSD = r_PSD_total(i)
              j = i
              
           end if
           
           
        end do
        
        
        buf(jpa,:,lpa,1) = Strauhal(j)
        buf(jpa,:,lpa,2) = maxPSD

	write(*,*)"FFT complete"


	Status = DftiFreeDescriptor(hand)


  






		
	
end subroutine
