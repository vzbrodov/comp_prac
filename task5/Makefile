FC = gfortran
FFLAGS = -O3

OBJ = main.o tridiagonal.o spline.o five_point_run.o

all: main

main: $(OBJ)
	$(FC) $(FFLAGS) -o main $(OBJ)

main.o: main.f95 spline.mod
	$(FC) $(FFLAGS) -c main.f95

tridiagonal.mod tridiagonal.o: tridiagonal.f95
	$(FC) $(FFLAGS) -c tridiagonal.f95

spline.mod spline.o: spline.f95 tridiagonal.mod five_point_run.mod
	$(FC) $(FFLAGS) -c spline.f95

five_point_run.mod five_point_run.o: five_point_run.f95
	$(FC) $(FFLAGS) -c five_point_run.f95

%.mod: %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod main

