comp=gfortran
main :main.o
	$(comp) main.f95 -o main -fopenmp
main.o:main.f95
	$(comp) -c $<
clean :
	rm *.o main
result : main
	./main
