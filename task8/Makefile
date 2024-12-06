FC = gfortran
FFLAGS = -O3

OBJ = main.o myfunction.o gauss_integr.o linear_system.o leg_polynomials.o

all: main

main: $(OBJ)
	$(FC) $(FFLAGS) -o main $(OBJ)

main.o: main.f95 gauss_integr.mod myfunction.mod
	$(FC) $(FFLAGS) -c main.f95

gauss_integr.mod gauss_integr.o: gauss_integr.f95 leg_polynomials.mod linear_system.mod
	$(FC) $(FFLAGS) -c gauss_integr.f95

%.mod: %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod main

