program one 
        use my_math
        implicit none

        integer, parameter :: q = 100
        integer, parameter :: IN = 1, OUT = 2
        real(8), parameter :: PI = 4.d0*datan(1.d0)
        integer :: N, i, Ns
        real(8) :: a, b, h
        real(8), dimension (:,:), allocatable :: grid
        character(5) :: method 

        call getarg (1, method) 

        select case (method)

        case default

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
                        grid(1, i) = a + i*h
                        grid(2, i) = cheb_interp(grid(1, i))
                enddo 

                open(OUT, file = "out/res_chebyshev.dat")  
                        write (OUT, "(2(f16.8, 2x))") grid(:,:)
                close(OUT)

        case ("uni")

                open (IN, file = 'in/uniform.dat')  
                        read (IN, '(2x, i10)') N
                        read (IN, *) a, b
                close (IN)

                Ns=N*q

                allocate (grid(1:2, 0:Ns))     

                h = (b-a)/ Ns

                do i = 0, Ns                           
                        grid(1, i) = a + i*h
                        grid(2, i) = unif_interp(grid(1, i))
                enddo 

                open(OUT, file = "out/res_uniform.dat")
                        write (OUT, "(2(f16.8, 2x))") grid(:,:)
                close(OUT)

        case ("cheb")

                open (IN, file = 'in/chebyshev.dat')       
                        read (IN, '(2x, i10)') N
                        read (IN, *) a, b
                close (IN)

                Ns=N*q

                allocate (grid(1:2, 0:Ns))   

                h = (b-a)/ Ns 

                do i = 0, Ns          
                        grid(1, i) = a + i*h
                        grid(2, i) = cheb_interp(grid(1, i))
                enddo 

                open(OUT, file = "out/res_chebyshev.dat")  
                        write (OUT, "(2(f16.8, 2x))") grid(:,:)
                close(OUT)

        end select

        deallocate(grid)


end program
