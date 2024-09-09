FC = gfortran
FFLAGS = -O3

OBJ = main.o function.o linear_system.o nonlinear_system.o

all: main

main: $(OBJ)
	$(FC) $(FFLAGS) -o main $(OBJ)

main.o: main.f95 function.mod nonlinear_system.mod
	$(FC) $(FFLAGS) -c main.f95

linear_system.mod linear_system.o: linear_system.f95
	$(FC) $(FFLAGS) -c linear_system.f95

nonlinear_system.mod nonlinear_system.o: nonlinear_system.f95 linear_system.mod
	$(FC) $(FFLAGS) -c nonlinear_system.f95

function.mod function.o: function.f95
	$(FC) $(FFLAGS) -c function.f95

%.mod: %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod main

