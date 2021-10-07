subroutine readfile

use variables
!$ use omp_lib
implicit none

real(8) :: c


!------------program start

!----read grid

!open(10,file='grid/grid.000001',form='unformatted')
!read(10)jmax,kmax,lmax

!write(*,*)'grid point',jmax,kmax,lmax


!allocate(x(jmax,kmax,lmax))
!allocate(y(jmax,kmax,lmax))
!allocate(z(jmax,kmax,lmax))

!write(*,*)'ndim	mdim'
!write(*,*)ndim,mdim

!read(10)x,y,z

close(10)

!write(*,*)jmax,kmax,lmax
!----------------------------
!-------read input
	a0 = 343.0d0
	theta = 0.0d0
	uj = a0 * 0.9d0
	rho0 = 1.0d0

	pi = 3.14159265358979323846d0
	theta = theta * 2.0d0*pi/360.0d0
	write(*,*) cos(theta),sin(theta)
   
    ! eg for Tcos0 (change for Tcos1 ...)
    ! Instantaneous fields

  write(filename,'("./CBogey2020/Jetv01D0420Re3125BL/signinstT0light.bin")')
  open(804,file=filename,form='unformatted')
  write(*,*)"open",filename
   ! open(804,file='./CBogey2020/Jetv01D0420Re3125BL/signinstTcos0light.bin',form='unformatted',status='unknown')
!	write(*,*)"open ./CBogey2020/Jetv01D0420Re3125BL/signinstTcos0light.bin"
    rewind(804)
    read(804) nbinst                ! nb of sampings
    read(804) nr                    ! nb of points in the radial direction
    read(804) nz                    ! nb of points in the axial direction
    read(804) thetagT0              ! azimuth
    
    print *, nbinst,nr,nz
    allocate(rg(nr))
    allocate(zg(nz))
    read(804) (rg(i),i=1,nr)        ! grid in the radial direction
    read(804) (zg(n),n=1,nz)        ! grid in the axial direction
    read(804) deltat
    write(*,*) deltat                ! time step
    read(804) r0                    ! jet radius
    write(*,*) r0
    deltat_nondim = deltat/(2.0*r0/uj)
    write(*,*) 'deltat_nonidim=',deltat_nondim 
    
   ! print *, rg/r0
   ! print *, zg/r0
 
    nb = 1 
    write(*,*)'flow import start'
    
    !------
	jmax = nz
	kmax = nb
	lmax = nr
	!------


	allocate(x(jmax,kmax,lmax))
	allocate(y(jmax,kmax,lmax))
	allocate(z(jmax,kmax,lmax))

	!data for cylindrical coordinate (instantaneous)

	 allocate(rhoT0(nz,nr))
    allocate(urT0(nz,nr))
    allocate(utT0(nz,nr))
    allocate(uzT0(nz,nr))
    allocate(pT0(nz,nr))
	allocate(vortT0(nz,nr))
	allocate(dilatT0(nz,nr))

	!-----------------------------------------------

	!-------------------------------------------------
	!data for cartecian cordinate (instantaneous)
	allocate(rho(jmax,lmax,nbinst))
	allocate(u(jmax,lmax,nbinst))
	allocate(v(jmax,lmax,nbinst))
	allocate(w(jmax,lmax,nbinst))
	!allocate(e(jmax,lmax))
	allocate(p(jmax,lmax,nbinst))
	!-------------------------------------------------


    ds = 1.0d0
    do j = 1,jmax      ! samplings
       do k = 1,kmax    ! axial direction
          do l = 1,lmax ! radial direction
             x(j,k,l) = zg(j)/r0
             y(j,k,l) = rg(l) * dsin(0.d0) + ds*dble(ii)
             z(j,k,l) = rg(l) * dcos(0.d0)/r0 
          enddo
       enddo
    enddo
	interval = 1
	
    open(80,file="case0.g",form="unformatted")
    write(80)jmax,kmax,lmax
        write(80) (((x(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),&
         (((y(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax),&
         (((z(j,k,l),j=1,jmax),k=1,kmax),l=1,lmax)
    close(80)





	dx = x(2,1,1) - x(1,1,1)
	dz = z(1,1,2) - z(1,1,1)

   do ii=1,nbinst !nbinst
	 nbg = ii*interval
    write(*,*)nbg,'/',nbinst,'ok'
    ! density, velocity components and pressure, vorticity norm and pressure

       read(804) ((rhoT0(j,l),l=1,nr),j=1,nz)
       read(804) ((urT0(j,l),l=1,nr),j=1,nz)
       read(804) ((utT0(j,l),l=1,nr),j=1,nz)
       read(804) ((uzT0(j,l),l=1,nr),j=1,nz)
       read(804) ((pT0(j,l),l=1,nr),j=1,nz)
       !read(804) ((vortT0(j,l),l=1,nr),j=1,nz) !for Tcos0
       !read(804) ((dilatT0(j,l),l=1,nr),j=1,nz)!for Tcos0

     !!! here your may need to store the fields in 3D tables (nr,nz,nbtinst)
      
       do j=1,nz
       do l=1,nr
    
    !--------------Nondimensionalization------------
     
          rhoT0(j,l) = rhoT0(j,l) / rho0
          urT0(j,l) = urT0(j,l) / uj
          utT0(j,l) = utT0(j,l) / uj
          uzT0(j,l) = uzT0(j,l) / uj
          pT0(j,l) = pT0(j,l) / (rho0 * uj * uj)
    
    !-----------------------------------------------

    !------cartesian verocity u v w 
          uu = (uzT0(j,l))
          vv = ( rg(l)*utT0(j,l)*dcos(0.0d0) + urT0(j,l)*dsin(0.0d0) )
          ww = ( urT0(j,l)*dcos(0.0d0) -rg(l)*utT0(j,l)*dsin(0.0d0) )
		   !vv = 0.0e0
    !!------------------------

          rho(j,l,ii) = rhoT0(j,l) 
          u(j,l,ii) = uu  
          v(j,l,ii) = vv
          w(j,l,ii) = ww
          p(j,l,ii) = pT0(j,l)
          

       end do
    	end do

	end do

	

end subroutine


















