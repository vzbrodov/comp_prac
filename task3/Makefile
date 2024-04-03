comp=gfortran
main: main.o linear_system.o
	$(comp) -O $^ -o $@ -fcheck=mem
main.o: main.f95 linear_system.mod
	$(comp) -O -c $<
linear_system.mod linear_system.o: linear_system.f95
	$(comp) -O -c linear_system.f95
clean:
	rm -f *.o main *.mod 	
result: main
	./main
res_gauss: main 
	./main g
res_jordan: main 
	./main j
res_choice: main 
	./main c
