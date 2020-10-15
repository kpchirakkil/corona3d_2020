corona3d_2020 aims to be a flexible planetary exosphere modeling program, able to handle
multiple particle types with a variety of initial distribution options and background
species profiles. See the sample 'corona3d_2020.cfg' file under the 'src' folder for
details on the currently available configuration parameter options.

This software is still in the early development stages, and there is much work yet to
be done. Below is a non-comprehensive 'to-do' list:

- change particle deactivation to particle replacement based on still-to-be-determined criteria
- continue developing Hot_H distribution class (currently includes H + H+ -> H+ + H* and HCO+ DR)
- figure out how to calculate and record escape flux of tracked particles
- consider adding inelastic collisions for molecular collision partners (all collisions are currently elastic)