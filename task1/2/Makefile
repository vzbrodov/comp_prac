comp=gfortran
main: main.o tridiagonal.o
	$(comp) $^ -o $@
main.o: main.f95 tridiagonal.mod
	$(comp) -c $<
tridiagonal.mod tridiagonal.o: tridiagonal.f95
	$(comp) -c tridiagonal.f95
clean:
	rm -f *.o main *.mod 	
result: main
	./main
