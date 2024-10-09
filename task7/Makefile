FC = gfortran
FFLAGS = -O3

OBJ = main.o fast_fourier.o

all: main

main: $(OBJ)
	$(FC) $(FFLAGS) -o main $(OBJ)

main.o: main.f95 fast_fourier.mod
	$(FC) $(FFLAGS) -c main.f95

fast_fourier.mod fast_fourier.o: fast_fourier.f95
	$(FC) $(FFLAGS) -c fast_fourier.f95

%.mod: %.f95
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod main

