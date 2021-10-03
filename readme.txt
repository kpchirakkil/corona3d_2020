corona3d_2020 aims to be a flexible planetary exosphere modeling program, able to handle
multiple particle types with a variety of initial distribution options and background
species profiles. See the sample 'corona3d_2020.cfg' file under the 'src' folder for
details on the currently available configuration parameter options.

**********************
** QUICKSTART GUIDE **
**********************

This software has been tested on both Linux and Mac systems using the GCC compilier.
To get started, first download the source code to the desired directory from a terminal:

     cd <desired directory>
     git clone https://github.com/rdelliott76/corona3d_2020.git

Now cd into the src directory and compile the application:

     cd corona3d_2020/src/
     make

Note that you may receive a compilation error due to missing the Eigen library. If this
happens, try this on Ubuntu:

     sudo apt install libeigen3-dev

On Mac:

     brew install eigen

You may need to look up how to install the Eigen library on your particular system if
neither of these commands works. Once installed, retry make, and it should work now.

Once installed, modify the corona3d_2020.cfg file to have the desired parameters. To
run the program, execute the following command from within the src directory:

     ./corona3d_2020



*************************
** SOURCE CODE DETAILS **
*************************

Following is a (somewhat) detailed explanation of the source code files and program logic
to aid future developers who would like to contribute to this software. There are also
fairly extensive comments throughout all of the source code files to aid in understanding.
A general knowledge of object-oriented programming concepts in C++ will be immensely helpful.

1. General Program Logic Overview

The first noteable bit of code that gets executed upon running the program is in the
Common_Functions.cpp file, specifically the get_seed() function. This bit of code seeds the
random number generator. It first looks for a file named "rng_seed" within the same directory
as the executable. This file should contain a single long long integer if you would like to run
a replicateable simulation for testing purposes. If the "rng_seed" file does not exist, then
the current system clock time is used as the seed. Execution now passes to main.cpp.

The main.cpp source file begins by reading in the main configuration file. This file must be
named "corona3d_2020.cfg" and must be in the same directory as the executeable.

*remove avg_thermal_v from main and atmosphere




