program main
    use :: myfunction
    use :: nonlinear_system
implicit none

    integer, parameter :: IN = 1, OUT = 2, maxiter = 100
    integer :: n, i
    real, allocatable :: x(:),res(:)

    open (IN, file = 'in/data.dat')        
    read (IN, '(2x, i10)') n ! размерность вектора X, пишется также в начале файла перед решеткой #
    allocate (x(1:n))

    do i=1,n
        read(IN, *) x(i)
    enddo
    close (IN)  

    res = newton(f, x, maxiter)
    write(*,*)"Модуль вектора F(X) = ", dot_product(f(res), f(res))

    open (OUT, file = "out/result.dat")
        write(OUT, "(f16.8)") res
    close(OUT)

deallocate(res, x)

end program main
