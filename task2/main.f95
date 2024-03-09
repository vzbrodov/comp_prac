program one 
        use my_math
        implicit none

        integer, parameter :: q = 100
        integer, parameter :: IN = 1, OUT = 2
        real, parameter :: PI = 4.d0*datan(1.d0)
        integer :: N, i, Ns
        real :: a, b, h
        real, dimension (:,:), allocatable :: grid 


        open (IN, file = 'in/uniform.dat')  
                read (IN, '(2x, i10)') N
                read (IN, *) a, b
                close (IN)

                Ns=N*q

                allocate (grid(1:2, 0:Ns)) 

                h = (b-a)/ Ns

                do i = 0, Ns                   ! равномерная сетка          
                        grid(1, i) = a + i*h
                        grid(2, i) = unif_interp(grid(1, i))
                enddo 

                open(OUT, file = "out/res_uniform.dat")
                        write (OUT, "(2(f16.8, 2x))") grid(:,:)
                close(OUT)

                do i = 0, Ns                    ! чебышевская сетка
                        grid(1, i) = ((a+b) + (b-a)*cos( (2.0*i + 1.0)*pi/(2.0*n + 2.0) ))/2.0
                        grid(2, i) = cheb_interp(grid(1, i))
                enddo 

                open(OUT, file = "out/res_chebyshev.dat")  ! вывод дискретизации графика сетки Чебышева в файл
                        write (OUT, "(2(f16.8, 2x))") grid(:,:)
                close(OUT)


        deallocate(grid)


end program
