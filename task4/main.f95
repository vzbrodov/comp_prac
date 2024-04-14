program main
    use :: iterative_methods_sle
implicit none

    integer :: n
    integer, parameter :: IN = 1, OUT = 2
    character(5) :: method !j-jacobi, z-zeidel, r - relaxation
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
        case default 
            write(*, *)"Неправильный ключ. Нужно выбрать: 'j'-jacobi, 'z'-zeidel, 'r' - relaxation"
            stop 0 
        case("j")
            X = jacobi (A, B)
            write(*, *)"Метод Якоби"
        case("z")
            X = zeidel (A, B)
            write(*, *)"Метод Зейделя"
        case("r")
            X = relaxation (A, B)
            write(*, *) "Метод релаксации"
        end select

        discr = sqrt(sum((matmul(A, X) - B)**2))
        write(*, *)'Модуль вектора невязок:', discr
            
        open (OUT, file = "out/result.dat")
            write(OUT,  "('# ', i10)") n
            write(OUT,  "(f16.7)" ) X
        close (OUT)

deallocate (A, B, X)

end program main
