module spline
    use :: tridiagonal
    use :: five_point_run
implicit none

contains

subroutine create_QAB (X, P, A, B, Q)
    implicit none
    real :: X(0:), P(0:)
    real :: A(1:,0:), B(1:,0:), Q(1:,0:) 
    integer :: n, i
    real :: l, c, r ! чтобы по одному разу считать элементы для матриц

    n = size(X) - 1

    A = 0.0
    B = 0.0
    A(2, 0) = 2*(X(1)-X(0))

    do i = 1, n-1
        l = X(i) - X(i-1)
        c = X(i+1) - X(i-1)
        r = X(i+1) - X(i)

        A(1,i) =  l
        A(2,i) = 2*c
        A(3,i) = r

        B(1, i) = 1/l
        B(2, i) = -(1/l)-(1/r)
        B(3, i) = 1/r
    enddo

    A(1, 1) = 0 
    A(3,n-1) = 0
     A(1, n) = 0
        A(3, n) = 0
    A(2, n) = 2*(X(n) - X(n-1))

    Q = 0
        do i = 0, n
            Q(2, i) = 1/P(i)
        enddo 


end subroutine create_QAB

subroutine create_SR (A, B, Y, Q, S, R)
    implicit none
    real :: A(1:,0:), B(1:,0:), Q(1:,0:)
    real ::Y(0:), S(0:), R(0:)
    real :: Bt (3,0:(size(S)-1)), matrix_1(5,0:(size(S)-1)) !Bt - транспонированная к B
    real :: matrix_2(0:(size(S)-1)) !matrix_2=6*B*Y
    integer :: i, n               !matrix_1=(A+6*B*Q*Bt) 

    n = size(S) - 1

        ! ищем транспонированную матрицу к B


    Bt(2,0) = B(2,0)
    Bt(3, 0) = B(1, 1)

    do i = 1, n
        Bt(1, i) = B(3,i-1) 
        Bt(2, i) = B(2, i)
        Bt(3, i) = B(1,i+1)
    enddo 

    matrix_1 = tridiag(B,Q) 
    matrix_1 = tridiag(matrix_1(2:4,:),Bt)
    matrix_1 = matrix_1*6
    matrix_1(2:4,:) = matrix_1(2:4,:) + A

    matrix_2 = 0
    do i = 1, n-1
        matrix_2(i) = 6*dot_product(B(:, i),Y(i-1:i+1))
    enddo

    S = fivepoint_run_method (matrix_1, matrix_2) 

    matrix_1 = tridiag(Q,Bt) !matrix_2=Y - Q*Bt*S
                              !matrix_1=Q*Bt  
    matrix_2 = 0

    matrix_2(0) = dot_product(matrix_1(3:5,0),S(0:2))
    matrix_2(1) = dot_product(matrix_1(2:5,1),S(0:3))

    do i = 2, n-2
        matrix_2(i) = dot_product(matrix_1(:,i),S(i-2:i+2))        
    enddo

    matrix_2(n-1) = dot_product(matrix_1(1:4,n-1),S(n-3:n))
    matrix_2(n) = dot_product(matrix_1(1:3,n),S(n-2:n))
    R = Y - matrix_2
end subroutine create_SR

function spline_f (chi, X, S, R) result (res)
    real :: chi, res
    real :: S(0:), R(0:), X(0:)
    integer ::  i, n
    real :: t, h
        
    n = size(X)-1
    i = minloc(abs(X-chi), dim = 1) - 1
    if (X(i) >= chi .and. i > 1) then 
        i = i-1
    endif
    h = X(i+1) - X(i)
    t = (chi - X(i))/h
    res = R(i)*(1-t) + R(i+1)*t - h*h*t*((1-t)/6 )*((2-t)*S(i) + (1+t)*S(i+1))
end function spline_f

end module spline
