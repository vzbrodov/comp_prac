comp=gfortran
main: main.o my_math.o
	$(comp) $^ -o $@
main.o: main.f95 my_math.mod
	$(comp) -c $<
my_math.mod my_math.o: my_math.f95
	$(comp) -c my_math.f95
clean:
	rm -f *.o main *.mod 	
result: main
	./main
