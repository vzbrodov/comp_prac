comp=gfortran
main: main.o iterative_methods_sle.o
	$(comp) -O $^ -o $@ -fcheck=mem
main.o: main.f95 iterative_methods_sle.mod
	$(comp) -O -c $<
iterative_methods_sle.mod iterative_methods_sle.o: iterative_methods_sle.f95
	$(comp) -O -c iterative_methods_sle.f95
clean:
	rm -f *.o main *.mod 	
result: main
	./main
jacobi: main 
	./main j
zeidel: main 
	./main z
relaxation: main 
	./main r
