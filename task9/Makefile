FC = gfortran
FFLAGS = -O3

OBJ = main.o settings.o diff_eq.o gauss_integr.o linear_system.o leg_polynomials.o nonlinear_system.o

all: main

main: $(OBJ)
	$(FC) $(FFLAGS) -o main $(OBJ)

main.o: main.f95 diff_eq.o settings.o
	$(FC) $(FFLAGS) -c main.f95

diff_eq.mod diff_eq.o: diff_eq.f95 gauss_integr.mod settings.mod nonlinear_system.mod
	$(FC) $(FFLAGS) -c diff_eq.f95

gauss_integr.mod gauss_integr.o: gauss_integr.f95 leg_polynomials.mod linear_system.mod 
	$(FC) $(FFLAGS) -c gauss_integr.f95

nonlinear_system.mod nonlinear_system.o: nonlinear_system.f95 linear_system.o
	$(FC) $(FFLAGS) -c nonlinear_system.f95

%.mod: %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod main

