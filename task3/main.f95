program main
    use linear_system
    implicit none

    integer :: n
    integer, parameter :: IN = 1, OUT = 2
    character(5) :: method !j-jordan, g-gauss, c - choice the leading element
    real, allocatable :: A(:,:), B(:), X(:)
    real :: discr !discrepancies

    call getarg (1, method) 

    open (IN, file = 'in/data.dat')        
    read (IN, '(2x, i10)') n 
    allocate (A(n, n), B(n), X(n))
    read (IN, *) A 
    read (IN, *) B 
    close (IN)  
    A = transpose(A) 
    select case (method)
        case ("g") !The basic scheme of the Gauss method
            X = gauss (A, B)
        case ("j") ! Jordan's scheme
            X = jordan (A, B)
        case ("c") !the scheme with the choice of the leading element
            X = choice (A, B)
    end select

    open (IN, file = 'in/data.dat')        
    read (IN, '(2x, i10)') n 
    read (IN, *) A 
    read (IN, *) B 
    close (IN)  
    discr = sqrt(sum((matmul(A, X) - B)**2))

    write(*, *)'Модуль вектора невязок:', discr

        open (OUT, file = "out/result.dat")
            write(OUT,  "('# ', i10)") n
            write(OUT,  "(f16.2)" ) X
        close (OUT)

    deallocate (A, B, X)
end program
