pure function bg_density(n_ref,h,r)
	use planet
	double precision bg_density
	double precision, intent(in) :: n_ref,h,r
	bg_density = n_ref*exp((r_ref-r)/h)
	return
end function 

program main

	use const
	use planet
	use phys

	implicit none

	integer*4 :: i, j, k, igen, ios, ix, jx, kx, ir, nv, nr, imain
	integer*8, dimension(-256:256,-256:256) :: density
	double precision :: ev(3), vtarg(3), vsec
	double precision, dimension(256,-90:90) :: poldens
	double precision, dimension(8) :: column
	double precision, dimension(-10:10,2) :: mu_r
	double precision :: samp, dv, u, v, rvi, mu, v_r, v_r_avg, v_avg, v_par, v_par_avg, g_asym, tau, u1(3), u2(3), m1, m2
	double precision :: rslop, r_exo, Odens, COdens, N2dens, CO2dens, mtarget, dens_conv, v_Obg
	parameter (nr = 49999)
	parameter (nv=69)
	parameter (dv=0.1)
	integer*8 vbin(nv),nvb
	integer*8, dimension(nr) :: obs_col, obs_dens
	character(len=1) :: istr

	double precision, external :: bg_density
	external collider

!	*** Initialization ***

	write(*,*) 'Initializing Simulation…'

	call random_seed()

	open(11,file='sim_param.list',iostat=ios,action='read')
	if (ios .ne. 0) then
		write(*,*) 'WARNING: Simulator parameter file is absent or bad. Using default values.'
	else
		read(11,nml=sim_param)
	endif
	close(11)

	open(11,file='fit_param.list',iostat=ios,action='read')
	if (ios .ne. 0) then
		write(*,*) 'WARNING: Fitting parameter file is absent or bad. Using default values.'
	else
 		read(11,nml=fit_param)
	endif
	close(11)

	density(:,:) = 0
	poldens(:,:) = 0
	column(:) = 0
	obs_col(:) = 0
	obs_dens(:) = 0
	mu_r(:,:) = 0
	samp = 0
	vbin = 0
	v_r = 0
	v_r_avg = 0
	v_avg = 0
	v_par = 0
 	v_par_avg = 0

	v_Obg = sqrt(8d0*k_b*T_bg/(pi*mbgsp(1)))	! velocity of background thermal O, used for calculating if hot test particles should be deactivated

	rslop = G*m_planet*dt**2/(4*r_b**2)		! slop in rentry radius due to timestep discretization, to prevent grazing eternal orbits

	call load_diff_xsect()
	call init_planet()
	call init_parts()
	
	! output distribution of positions and velocities to file
    open(unit = 55, file='positions.out')
    open(unit = 56, file='velocities.out')
    do j = 1 , n_part_max
       write(55,*) parts(j)%x(1), parts(j)%x(2), parts(j)%x(3)
       write(56,*) parts(j)%v(1), parts(j)%v(2), parts(j)%v(3)
    enddo
    close(unit=55)
    close(unit=56)

  write(*,*) 'Simulating Particle Transport...'
 
!	*** Main Loop  ***
  do imain = 1,nstep!nstep_eq
!  do while (n_escape .lt. nstep_err)
!  	write(*,*) imain
	if ( mod(imain,10000) .eq. 0) write(*,'(2I16)') imain, sum(parts(:)%active)
	call step(dt)

	do j=1,n_part

!	*** Cycle if ith Test Particle is Deactivated ***
		if (parts(j)%active .eq. 0) cycle

!	*** Calculate Magnitude of Test Particle Velocity Vector ***
		v = sqrt(sum(parts(j)%v(:)**2))

!	*** Record in Integrated Column Density ***
		ix = int( (parts(j)%x(1) - r_planet)/1d3 )
		if ( ((ix.ge.300).and.(ix.le.nr)) .and. (abs(parts(j)%x(2)).le.500d3)  ) then
!			obs_col( ix ) = obs_col( ix ) + parts(j)%w	! if you are weighting the particles, you would record them this way
			obs_col( ix ) = obs_col( ix ) + 1
		endif

!	*** Record in Day Hemisphere Average Radial Density ***
		ix = int( (parts(j)%r - r_planet)/1d3 )
                if ( ((ix.ge.300).and.(ix.le.nr)) .and. (parts(j)%x(1).ge.0d0)  ) then
                        obs_dens( ix ) = obs_dens( ix ) + 1
                endif

		! Assume symmetry for a given solar zenith angle
                ix = nint(parts(j)%x(2)/200d3)
                jx = nint(parts(j)%x(1)/200d3)
		kx = nint(parts(j)%x(3)/200d3)
                if ((abs(kx).le.256).and.((abs(ix).le.256).and.(abs(jx).le.256))) then
			density( ix, jx ) = density( ix, jx ) + 1
			density(-ix, jx ) = density(-ix, jx ) + 1
			density( kx, jx ) = density( kx, jx ) + 1
			density(-kx, jx ) = density(-kx, jx ) + 1
                endif

!	*** Example for Recording Velocity Distribution ***
!		if ( ( abs(parts(j)%r - (r_planet+500d3)) .lt. 10d3) .and. (parts(j)%x(1) .gt. 0d0) ) then		! record particles on dayside within 10 km of 500 km altitude 
!!			write(*,*) v
!			v_r = dot_product(parts(j)%v,parts(j)%x)/parts(j)%r
!			v_avg = v_avg + v
!			v_r_avg = v_r_avg + abs( v_r )
!			v_par = sqrt((v**2-v_r**2)/2d0)
!			v_par_avg = v_par_avg + v_par
!
!			nvb = 1+int(v/1d3/dv)
!			vbin(nvb) = vbin(nvb) + 1
!!			ix = int( atan(parts(j)%v(1)/sqrt(parts(j)%v(2)**2+parts(j)%v(3)**2)) * 2/pi * 180.0 / 10.0 ) + 1
!			ix = int( (dot_product(parts(j)%v,parts(j)%x)/(v*parts(j)%r) + sign(0.1,v_r))*10 )
!			mu_r( ix, 1 ) = mu_r( ix, 1 ) + v
!			mu_r( ix, 2 ) = mu_r( ix, 2 ) + 1
!!			write(*,*) mu_r(ix,2)
!		endif

!		/// Simulate Collisions \\\

!		*** Sample Background Species Densities at Position of Particle ***
		Odens = bg_density(n_bg200(idxO),H_bg(idxO),parts(j)%r)
		COdens = bg_density(n_bg200(idxCO),H_bg(idxCO),parts(j)%r)
		N2dens = bg_density(n_bg200(idxN2),H_bg(idxN2),parts(j)%r)
		CO2dens = bg_density(n_bg200(idxCO2),H_bg(idxCO2),parts(j)%r)

! 		*** Determine if Test Particle Collided ***
		tau = v*dt*( O_sigma*Odens + CO_sigma*COdens + N2_sigma*N2dens + CO2_sigma*CO2dens )
		call random_number(u)
		if ( u .gt. exp(-tau) ) then

			n_collide = n_collide + 1

!			*** Pick Target Species Mass ***
			call random_number(u)
			if ( u .lt. Odens/(Odens+COdens+N2dens+CO2dens) ) then
				mtarget = mbgsp(idxO)			! collision with O
			else if ( u .lt. (Odens+N2dens)/(Odens+COdens+N2dens+CO2dens) ) then
				mtarget = mbgsp(idxN2)			! collision with N2
			else if ( u .lt. (Odens+COdens+N2dens)/(Odens+COdens+N2dens+CO2dens) ) then
				mtarget = mbgsp(idxCO)			! collision with CO
			else
				mtarget = mbgsp(idxCO2)			! collision with CO2
			endif

!			*** Begin Inverse CDF Sampling of Differential Cross-Section ***
			call random_number(u)
			k=1
			do while (cdf(1,k) .lt. u)
				k=k+1
			enddo
			mu = cdf(2,k)
!			*** End Inverse CDF Sampling ***

!			*** Sample 3D Velocity Vector for Collision Partner using Maxwellian at Temperature = T_bg ***
			call gen_mb(sqrt(k_b*T_bg/mtarget), vtarg)
!			vtarg = (/0d0,0d0,0d0/) ! test line for zero temperature background gas

!			*** Collide the Test Particle with the Collision Partner ***
			call collider(parts(j)%v,vtarg,mbgsp(idxO),mtarget,mu)

!			vsec = sqrt(sum(vtarg(:)**2)) ! the velocity of the collision partner can be used for other things

		endif

!		*** Deactivate Particles as Necessary ***
!		the thermal O velocity v_Obg is included to make sure we don't ignore the possibility that a later collision could boost the apex of the hot particle a little higher
!		this is a real effect that significantly impacts the profile of bound oxygen at lower altitudes
!		if ( (parts(j)%r .lt. (r_planet + 600d3)) .and. (v + v_Obg .lt. sqrt(2d0*G*m_planet*(parts(j)%ri-1d0/(r_planet+600d3)))) ) then	! deactivate if it can't get above 600 km
		if ( (parts(j)%r .lt. (r_planet + 900d3)) .and. (v + v_Obg .lt. sqrt(2d0*G*m_planet*(parts(j)%ri-1d0/(r_planet+900d3)))) ) then	! deactivate if it can't get above 900 km

			parts(j)%active = 0
		endif
	enddo
  enddo
  
  ! output final positions and velocities
  open(unit = 75, file='positions2.out')
  open(unit = 76, file='velocities2.out')
  do j = 1 , n_part_max
     if (parts(j)%active .eq. 0) cycle
     write(75,*) parts(j)%x(1), parts(j)%x(2), parts(j)%x(3)
     write(76,*) parts(j)%v(1), parts(j)%v(2), parts(j)%v(3)
  enddo
  close(unit=75)
  close(unit=76)

dens_conv = dt*2d0*H_DR*DR160*2d0*pi*(r_b+H_DR)**2        				! density conversion factor; 2 hot O per DR across day hemisphere

write(*,'(A21,I14)') '  Spawned Particles:',int(n_spawn,kind=8)
write(*,'(A21,I14)') '   Active Particles:',sum(parts(:)%active)
write(*,'(A21,I14)') 'Particle Collisions:',int(n_collide,kind=8)
write(*,'(A21,I14)') '     Actual Samples:',n_part*nstep
!write(*,'(A21,ES14.3)') '  Average v @ 500km:',v_avg/sum(vbin(:))
!write(*,'(A21,ES14.3)') '   Radial v @ 500km:',v_r_avg/sum(vbin(:))
!write(*,'(A21,ES14.3)') ' Parallel v @ 500km:',v_par_avg/sum(vbin(:))

write(*,*) 'Global DR Rate', H_DR*DR160*2d0*pi*(r_b+H_DR)**2

dens_conv = dt*2d0*H_DR*DR160*2d0*pi*(r_b+H_DR)**2        				! density conversion factor; 2 hot O per DR across day hemisphere

!open(unit = 44, file='density1d'//istr//'.out')						! for 'poor man's parallelization' you can tag the outputs of simultaneous runs in the output name
open(unit = 44, file='density1d.out')
do i=300,nr
       write(44,'(F10.1,2ES14.6)') dble(i), &						! the following need division by altitude range to get [m-3] & [m-2]
		obs_dens(i)/dble(n_spawn)*dens_conv/(2d0*pi*(r_planet+i*1d3)**2), &	! average dayside radial density 
		obs_col(i)/dble(n_spawn)*dens_conv/1000d3				! observation line-of-sight column integrated density; along slit integration
enddo
close(unit=44)

open(unit = 45, file='density2d.out')
do i=-256,256
	do j=-256,256
		write(45,'(I4,I5,I14)') i,j, density( i, j )
	enddo
enddo
close(unit=45)

end program main
