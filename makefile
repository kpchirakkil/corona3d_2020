FC = gfortran
#FLAGS = -O0 -traceback -check all #-pg -openmp -parallel -fpp -ipo
#FLAGS = -O2 -fpp -ipo #-mcmodel=large -shared-intel #-openmp -parallel

FLAGS = -O2

MC_OBJ = const.o planet2.o phys.o collider2.o main.o

mc : $(MC_OBJ)
	$(FC) $(MC_OBJ) $(FLAGS) -o mc

%.mod : %.o
	@if [! -f $@ ]; then \
		rm $< \
		$(MAKE) $< \
	fi

%.o : %.f*
	$(FC) $(FLAGS) -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.mod
	rm -f mc
