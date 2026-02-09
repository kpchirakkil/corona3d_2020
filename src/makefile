# Auto-detect OS and set compiler
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS - use clang++ with aggressive optimization
    CXX := clang++
    CFLAGS := -O3 -march=native -ffast-math -funroll-loops
else
    # Linux - use g++ with aggressive optimization
    CXX := g++
    CFLAGS := -O3 -march=native -ffast-math -funroll-loops -Wall
endif

corona3d_2020: Atmosphere.o Background_Species.o Common_Functions.o Distribution_Hot_H.o Distribution_Hot_O.o Distribution_Import.o Distribution_MB.o Distribution.o Interpolator.o main.o Particle_CO.o Particle_CO2.o Particle_H.o Particle_N2.o Particle_O.o Particle.o Planet.o vtally.o
	$(CXX) $(CFLAGS) Atmosphere.o Background_Species.o Common_Functions.o Distribution_Hot_H.o Distribution_Hot_O.o Distribution_Import.o Distribution_MB.o Distribution.o Interpolator.o main.o Particle_CO.o Particle_CO2.o Particle_H.o Particle_N2.o Particle_O.o Particle.o Planet.o vtally.o -o corona3d_2020

Atmosphere.o: Atmosphere.cpp
	$(CXX) $(CFLAGS) -c Atmosphere.cpp

Background_Species.o: Background_Species.cpp
	$(CXX) $(CFLAGS) -c Background_Species.cpp

# Compiled without -ffast-math so that NaN/Inf detection and creation work correctly.
# This file only runs at startup (CSV import); no performance impact.
CFLAGS_SAFE := $(filter-out -ffast-math,$(CFLAGS))
Common_Functions.o: Common_Functions.cpp
	$(CXX) $(CFLAGS_SAFE) -c Common_Functions.cpp

Distribution_Hot_H.o: Distribution_Hot_H.cpp
	$(CXX) $(CFLAGS) -c Distribution_Hot_H.cpp

Distribution_Hot_O.o: Distribution_Hot_O.cpp
	$(CXX) $(CFLAGS) -c Distribution_Hot_O.cpp

Distribution_Import.o: Distribution_Import.cpp
	$(CXX) $(CFLAGS) -c Distribution_Import.cpp

Distribution_MB.o: Distribution_MB.cpp
	$(CXX) $(CFLAGS) -c Distribution_MB.cpp

Distribution.o: Distribution.cpp
	$(CXX) $(CFLAGS) -c Distribution.cpp

Interpolator.o: Interpolator.cpp
	$(CXX) $(CFLAGS) -c Interpolator.cpp

main.o: main.cpp
	$(CXX) $(CFLAGS) -c main.cpp

Particle_CO.o: Particle_CO.cpp
	$(CXX) $(CFLAGS) -c Particle_CO.cpp

Particle_CO2.o: Particle_CO2.cpp
	$(CXX) $(CFLAGS) -c Particle_CO2.cpp

Particle_H.o: Particle_H.cpp
	$(CXX) $(CFLAGS) -c Particle_H.cpp

Particle_N2.o: Particle_N2.cpp
	$(CXX) $(CFLAGS) -c Particle_N2.cpp

Particle_O.o: Particle_O.cpp
	$(CXX) $(CFLAGS) -c Particle_O.cpp

Particle.o: Particle.cpp
	$(CXX) $(CFLAGS) -c Particle.cpp

Planet.o: Planet.cpp
	$(CXX) $(CFLAGS) -c Planet.cpp

vtally.o: vtally.cpp
	$(CXX) $(CFLAGS) -c vtally.cpp

clean:
	rm -f *.o
	rm -f corona3d_2020
