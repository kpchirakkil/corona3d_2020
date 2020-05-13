collider2.f90		subroutine that calculates final velocity vectors of two colliding particles given masses, initial velocity vectors, and a scattering angle
const.f90		module that contains physical constants; don't change this unless you are God (or his designated representative)
main.f90		contains the main loop that initializes and time steps the simulation, randomly colliding the test particles with background atmosphere
phys.f90		module that contains the particle type, generates the initial test particle velocity, and integrates the equations of motion
planet.f90		module that contains physical characteristics of the test particles and the planet's size/gravity/background atmosphere

sim_param.list		namelist that allows configuration of simulation parameters, e.g. timestep, number of steps
fit_param.list		namelist that allows configuration of ionospheric characteristics for generating test particles

KHARCHENKOMAVEN.TXT	input file containing differential scattering cross-section for O-O collisions
KHARCHENKOCDF.TXT	input file containing CDF of differential scattering cross-section for O-O collisions