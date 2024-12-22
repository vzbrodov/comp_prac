FC = gfortran
FFLAGS = -O3 -z execstack

OBJ = main.o my_prec.o settings.o diff_eq.o gauss_integr.o linear_system.o leg_polynomials.o nonlinear_system.o

all: main

main: $(OBJ)
	$(FC) $(FFLAGS) -o main $(OBJ)

main.o: main.f95 my_prec.o diff_eq.o settings.o
	$(FC) $(FFLAGS) -c main.f95

my_prec.o: my_prec.f95
	$(FC) $(FFLAGS) -c my_prec.f95

settings.o: settings.f95 my_prec.o
	$(FC) $(FFLAGS) -c settings.f95

diff_eq.o: diff_eq.f95 gauss_integr.o settings.o nonlinear_system.o linear_system.o my_prec.o  
	$(FC) $(FFLAGS) -c diff_eq.f95

gauss_integr.o: gauss_integr.f95 leg_polynomials.o linear_system.o my_prec.o  
	$(FC) $(FFLAGS) -c gauss_integr.f95

nonlinear_system.o: nonlinear_system.f95 linear_system.o my_prec.o  
	$(FC) $(FFLAGS) -c nonlinear_system.f95
	
linear_system.o: linear_system.f95 my_prec.o  
	$(FC) $(FFLAGS) -c linear_system.f95

leg_polynomials.o: leg_polynomials.f95 my_prec.o
	$(FC) $(FFLAGS) -c leg_polynomials.f95

clean:
	rm -f *.o *.mod main
