module phys

  use const
  use planet
  use omp_lib

  implicit none

  integer*4 :: n_part_max, n_part		! Number of simulated particles
  parameter ( n_part_max = 10000 )		! this should be large enough to allow good statistical sampling, but not so large it takes forever to run

! /// Define Particle Type \\\
  type :: particle
     integer			:: active	! flag for whether particle is active, i.e. should still be considered in the simulation 
!     double precision		:: w		! weighting, not utilized right now but can be handy
     double precision           :: r		! radius from center of planet			[m]
     double precision           :: ri		! inverse radius (for computational efficiency)	[m-1]	
     double precision           :: x(3)		! position vector				[m,m,m]
     double precision           :: v(3)		! velocity vector				[m/s,m/s,/ms]
  end type particle

! /// Make an Array of Particles \\\
  ! type(particle), private, allocatable       :: parts(:)	! sometimes it can be more memory efficient to allocate at run time 
  type(particle)                             :: parts(n_part_max)

  double precision :: Ekin, Epot, Etot
  integer*4        :: n_escape  = 0
  integer*4        :: n_spawn   = 0
  integer*4        :: n_reenter = 0
  integer*4        :: n_collide = 0

!  double precision :: n_escape = 0 ! integers make sense for counting particles, but use doubles when weighting particles
!  double precision :: n_spawn  = 0
!  double precision :: n_reenter= 0
!  double precision :: n_collide = 0

  double precision dxs(2,179)		! array to hold differential scattering cross-section
  double precision cdf(2,180)		! CDF of differential scattering cross-section, for sampling

contains

! /////   Subroutine to load differential cross-section (Kharchenko et al. 2000)   \\\\\
subroutine load_diff_xsect()
	open(1,file='KHARCHENKOMAVEN.TXT')
	read(1,*)
        read(1,*) dxs
        close(1)
	open(1,file='KHARCHENKOCDF.TXT')
	read(1,*) cdf
	close(1)
end subroutine load_diff_xsect

! /////   Subroutine to populate initial particle position/velocity distribution   \\\\\
subroutine init_parts()

    integer*4        :: j
    double precision :: randnum

    double precision :: s, ds
    double precision :: rb, rt, H, nb, nv, nvt

    ! populate initial particles
    n_part = n_part_max
    do j = 1, n_part_max

       call random_number(randnum)
       call new_particle(j)

       ! write(*,*) parts(i)%r, phi, theta
       ! write(*,*) parts(i)%x(1), parts(i)%x(i,2), parts(i)%x(i,3)
    enddo
    
  end subroutine init_parts

! /////   Subroutine to Create a New Particle   \\\\\
  subroutine new_particle(index)

    integer*4, intent(in) :: index
    integer*4 j,k, done
    double precision r

    double precision :: randnum, phi, u, v, rvi, mu, vimp, x, f, a, g, va, Ti, Te, m_ion, Ei
    double precision, dimension(3) :: ve, vion
    save vimp
    data vimp /0d0/

        Ti = 400d0
        Te = 1600d0
	m_ion = 31.9983d0*amu    ! O2+ ion mass

!	altitude distribution for O2+ dissociative recombination
	call random_number(randnum)
!	r = r_b - log(randnum)*H_DR
	r = r_planet + 160d3 - log(randnum)*H_DR

!      ### Global Source ###
       call random_number(randnum)
       phi = 2d0*pi*randnum
       call random_number(randnum)
       u = 2d0*randnum-1
       call random_number(randnum)
       parts(index)%r = r
       parts(index)%ri = 1/r
       parts(index)%x(1) = r*sqrt(1-u*u)*cos(phi)
       parts(index)%x(2) = r*sqrt(1-u*u)*sin(phi)
       parts(index)%x(3) = r*u
!      ### Hemispherical Adjustment For Dayside Photochemical Process ###
       if ( parts(index)%x(1) .lt. 0d0 ) parts(index)%x(1) = -parts(index)%x(1)

!	### point source test case ###
!	parts(index)%r = r
!	parts(index)%ri = 1/r
!	parts(index)%x(1) = r
!	parts(index)%x(2) = 0d0
!	parts(index)%x(3) = 0d0

!	### initialize to zero velocity ###
	parts(index)%v(:) = 0d0

!	### Select Hot O Electronic Channel ###
!	Probabilities and energies
!	Particles are activated per channel to allow testing the contribution of each population in the simulation

        call random_number(randnum)
        if (randnum .lt. 0.22d0) then
                Ei = 6.99d0*jev
!		parts(index)%w = 1d0 			! you could weight the particles from different channels here, if you wanted
		parts(index)%active = 1
	else if (randnum .lt. (0.22+0.42)) then
                Ei = 5.02d0*jev
		parts(index)%active = 1
        else if (randnum .lt. (0.22+0.42+0.31)) then
                Ei = 3.06d0*jev
		parts(index)%active = 1
        else
                Ei = 0.84d0*jev
		parts(index)%active = 1
        end if

!	### Thermal Rotational Energy of O2+ ###
        Ei = Ei + e_rot(3.35967d-23,Ti)

!	### Translational Energy per O Resulting from Dissociative Recombination of O2+ ###
        v = sqrt(Ei/m)

!	### spherically isotropic velocity vector ###
       call random_number(randnum)
       phi = 2*pi*randnum
       call random_number(randnum)
       u = 2*randnum-1
       call random_number(randnum)
       parts(index)%v(1) = v*sqrt(1-u*u)*cos(phi)
       parts(index)%v(2) = v*sqrt(1-u*u)*sin(phi)
       parts(index)%v(3) = v*u

!       ### Add initial ion and electron translational momentum ###
        va = sqrt(k_b*Ti/(m_ion))	! average thermal ion velocity
        call gen_mb(va,vion)		! sample from Maxwell-Boltzmann distribution
        va = sqrt(k_b*Te/m_e)		! average thermal electron velocity
        call gen_mb(va,ve)		! sample from Maxwell-Boltzmann distribution
        parts(index)%v(:) = parts(index)%v(:) + (m_ion*vion(:) + m_e*ve(:))/(m_ion+m_e) ! add it all up

!	n_spawn = n_spawn + parts(index)%w ! if you wanted to weight particles, you can assign it here
	n_spawn = n_spawn + 1

  end subroutine new_particle


!	/////	Iterate Equation of Motion using Velocity-Verlet Algorithm  \\\\\
  subroutine step(dt)

    integer*4        :: i
    double precision :: dt, a(3)

    do i=1,n_part
       if (parts(i)%active .eq. 0) cycle
       ! calculate acceleration at current position
       a(:) = k_g*parts(i)%x(:)*(parts(i)%ri)**3

       ! calculate next position and update particle
       parts(i)%x(:) = parts(i)%x(:) +  parts(i)%v(:)*dt + 0.5d0*a(:)*dt*dt
       parts(i)%r = sqrt(sum(parts(i)%x(:)**2))
       parts(i)%ri = 1/parts(i)%r

       ! calculate acceleration at next position
       a(:) = a(:) + k_g*parts(i)%x(:)*(parts(i)%ri)**3

       ! calculate next velocity using acceleration at next position and update particle
       parts(i)%v(:) = parts(i)%v(:) + 0.5d0*a(:)*dt
    enddo
  end subroutine step


! /// Sample from a Maxwell-Boltzmann Distribution \\\
  subroutine gen_mb(va, v)

  double precision, intent(in) :: va         ! va = sqrt(k_b*T/m), average velocity	[m/s]
  double precision, intent(out) :: v(3)      ! 3D Maxwell-Boltzman velocity		[m/s,m/s,m/s]

  double precision :: randnum1, randnum2, randnum3, randnum4

  call random_number(randnum1)
  call random_number(randnum2)
  call random_number(randnum3)
  call random_number(randnum4)

  randnum1 = va*sqrt(-2d0*log(1d0-randnum1))
  randnum2 = twopi*randnum2
  randnum3 = va*sqrt(-2d0*log(1d0-randnum3))
  randnum4 = twopi*randnum4

  v(1) = randnum1*cos(randnum2)
  v(2) = randnum1*sin(randnum2)
  v(3) = randnum3*cos(randnum4)

  end subroutine gen_mb

! /// Sample Rotational Energy from a Thermal Population of Rigid Rotators \\\
function e_rot(B,T)
        integer j, nj
        double precision, intent(in) :: B, T		! B = molecular moment of inertia, T = rotational temperature
        double precision e_rot, Qrot, Ptot, u, Pj
        parameter (nj = 200)                            ! highest rotational level considered

        Qrot = k_b*T/B + 1d0/3d0                        ! parition function using two terms of Euler-MacLaurin summation formula

        Ptot = 0d0
        call random_number(u)
        do j=0,nj
                Pj = (2d0*j+1d0)*exp(-j*(j+1)*B/(k_b*T))/Qrot
                Ptot = Ptot + Pj
                if (u .lt. Ptot) exit
        enddo

        e_rot = j*(j+1d0)*B
        return
end function e_rot

end module phys
