program main
    use :: spline
implicit none

    integer, parameter :: IN = 1, OUT = 2
    integer, parameter :: qq = 100 
    integer :: n, i
    real :: chi, h
    real, allocatable :: X(:), Y(:), P(:), R(:), S(:)
    real, allocatable :: A(:,:), B(:,:), approx_func(:,:), Q(:,:)

    open (IN, file = 'in/data.dat')        
    read (IN, '(2x, i10)') n 
    allocate (X(0:n), Y(0:n), P(0:n))
    do i=0,n
        read(IN, *)X(i), Y(i), P(i)
    enddo
    close (IN)  

    allocate(R(0:n),S(0:n))
    allocate(A(3, 0:n),B(3,0:n),Q(3, 0:n))
    
    call create_QAB(X,P,A,B, Q)
    call create_SR (A, B, Y, Q, S, R)

    deallocate (A, B, Q, Y, P)

    h = (X(n) - X(0))/ (qq*n)
    allocate  (approx_func(2, 0:qq*n))

    do i = 0, qq*n
        chi = X(0) + i*h
        approx_func(1, i) = chi
        approx_func(2, i) = spline_f(chi, X, S, R)
enddo 

    open (OUT, file = "out/result.dat")
    write(OUT, "(2(f16.8, 1x))") approx_func
    close(OUT)

deallocate(approx_func, X, R, S)

end program main
