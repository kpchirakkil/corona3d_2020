corona3d_2020 aims to be a flexible planetary exosphere modeling program, able to handle
multiple particle types with a variety of initial distribution options and background
species profiles. See the sample 'corona3d_2020.cfg' file under the 'src' folder for
details on the currently available configuration parameter options.

This software is still in the early development stages, and there is much work yet to
be done. Below is a non-comprehensive 'to-do' list:

- parameterize particle deactivation criteria (as of now it's hardcoded in the atmosphere class)
- expand and parameterize simulation output options
- continue developing Hot_H distribution class
- add ability to use custom density profiles for background species
- add ability to use atmosphere temperature profile
- possibly add ability to import hdf5 file types for above mentioned profiles
- consider adding inelastic collisions for molecular collision partners (all collisions are currently elastic)