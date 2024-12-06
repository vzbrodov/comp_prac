program main
    use :: gauss_integr
    use :: myfunction
implicit none
    integer :: n
    character(2) :: n1
    real :: a, b, result

    a = 3
    b = 6
    call getarg (1, n1)
    read(n1, "(i4)") n
    
    result = gauss_int(f,a,b,n)
    WRITE(*,*) 'RESULT =',result
end program main
