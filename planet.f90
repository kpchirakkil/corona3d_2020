module planet

  use const

  implicit none

  ! *** Simulation Settings ***
  double precision :: dt         	= 5.0d-1			! time step                     [s]
  double precision :: esc_err	 	= 0.01d0			! accuracy of escape rate
  integer*8        :: nstep      	= 100000			! number of steps in simulation

  double precision :: r_b		=  3557d3	       		! bottom of model		[m]
  double precision :: r_t		= 10000d3        	        ! top of model			[m]

  ! *** Planetary Parameters ***
  double precision, parameter :: r_planet	= 3397d3                ! planetary radius              [m]
  double precision, parameter :: m_planet	= 6.4185d23		! planetary mass                [kg]

  ! *** Thermospheric Species ***
  integer, parameter          :: nbgsp		= 4			! number of background species  
  double precision, parameter :: mbgsp(nbgsp)	= (/15.9994d0*amu, &	!   O mass			[amu]
						    28.0101d0*amu, &	!  CO mass
						    28.0134d0*amu, &	!  N2 mass
						    44.0095d0*amu  &	! CO2 mass
                                                	/)
!  double precision, parameter :: xsect(nbgsp)	= (/6.4d-19 &		! O cross-section               [m-2]
!                                                	/)
  double precision :: m         = 15.9994d0*amu				! escaping species mass         [kg]
  double precision :: O_sigma   = 6.4d-19				! O cross section               [m^2]
  double precision :: CO_sigma  = 1.85d-18				! CO cross section              [m^2]
  double precision :: N2_sigma  = 1.85d-18				! N2 cross section              [m^2]
  double precision :: CO2_sigma = 2.0d-18				! CO2 cross section             [m^2]

  ! *** Thermospheric Structure ***
  double precision, parameter :: r_ref		= r_planet + 200d3	! reference altitude		[m]
  double precision :: T_bg			= 277.6d0		! background temperature	[K]
  double precision :: H_bg(nbgsp)					! scale heights			[m]

  double precision :: n_bg200(nbgsp)		= (/2.64d13, &		! O density @ 200 km		[m-3]
						    1.00d13, &		! CO
						    3.10d13, &		! N2
						    6.68d13  &		! CO2
						/)

  double precision :: DR160		= 600d6				! Dissociation Recombination	[m-3/s]
									! 	rate at 200 km
  double precision :: H_DR		= 22d3				! DR scale height		[m]

  double precision k_g, vt0
  integer*8 nstep_err

  integer idxO, idxCO, idxN2, idxCO2
  parameter (idxO   = 1)
  parameter (idxCO  = 2)
  parameter (idxN2  = 3)
  parameter (idxCO2 = 4)

  namelist /sim_param/ dt, esc_err, nstep, nstep_err, r_b, r_t
  namelist /fit_param/ T_bg, n_bg200, DR160, H_DR

contains

! intialize the simulation parameters
  subroutine init_planet

    k_g = -G*m_planet

    H_bg(:) = k_b*T_bg/(mbgsp(:)*3.31d0)

    vt0 = sqrt( k_b*T_bg/m )

    nstep_err = nint(1/esc_err**2)

  end subroutine init_planet

end module planet
