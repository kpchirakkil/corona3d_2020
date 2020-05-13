module const

  implicit none

  double precision pi, twopi, k_b, c, G, amu, m_e, q_e, jev

  parameter (pi    = 3.1415926535897932d0)              ! pi                            unitless
  parameter (twopi = 2*pi)                              ! 2*pi                          unitless
  parameter (k_b   = 1.3806504d-23)                     ! Boltzmann's Constant          [J/K]
  parameter (c     = 299792458d0)                       ! Speed of Light in Vacuum      [m/s]
  parameter (G     = 6.67428d-11)                       ! Gravitational Constant        [m^3/kg/s^2]
  parameter (amu   = 1.660538782d-27)                   ! Atomic Mass Unit              [kg]
  parameter (m_e   = 9.10938215d-31)                    ! Electron Mass                 [kg]
  parameter (q_e   = 1.602176487d-19)                   ! Elementary Charge             [C]
  parameter (jev   = q_e)                               ! Joules/Electron Volt          unitless

end module const
