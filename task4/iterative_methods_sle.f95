module iterative_methods_sle
    implicit none
    real, parameter :: eps = 1e-6
contains

function jacobi (A, B) result (X)
    implicit none
    real :: A(1:, 1:), B(1:)
    real :: X(1:size(B))
    real :: D(1:size(B), 1:size(B)), inverseD(1:size(B), 1:size(B)), Z(1:size(B), 1:size(B)) !inverseD - обратная к D
    real :: Xk(1:size(B)), G(1:size(B))
    integer :: n, i, j
    
    do i = 1, n 
        if(abs(A(i, i)) < sum(abs(A(i,1:i-1))) + sum(abs(A(i,i+1:n)))) then       ! проверка на диагональное преобладаение
             write(*, "('Матрица не обладает диагональным преобладанием')")
        endif
    enddo

    n = size(B)
    D = 0
    inverseD = 0 

    do i = 1, n
        D(i,i) = A(i,i)
        inverseD(i,i) = 1/A(i,i)
    enddo 

    Z = matmul(inverseD, D - A)
    G = matmul(inverseD, B)

    X = 0
    Xk = G

    do while (sqrt(sum((Xk - X)**2)) >= eps)
        X = Xk
        Xk = matmul(Z, X) + G
    enddo

    X = Xk
end function jacobi



function zeidel (A, B) result (X)
    implicit none
    real :: A (1:, 1:), B (1:)
    real :: X (1:size(B))
    real :: P(1:size(B), 1:size(B))
    real :: Q(1:size(B)), Xk(1:size(B))
    integer :: n, i, j

    do i = 1, n 
        if(abs(A(i, i)) < sum(abs(A(i,1:i-1))) + sum(abs(A(i,i+1:n)))) then       ! проверка на диагональное преобладаение
            write(*, "('Матрица не обладает диагональным преобладанием')")
        endif
    enddo

    n = size (B)

    do i=1,n
        do j=1,n
            P(i,j) = - A(i,j)/A(i,i)
        enddo
    enddo

    do i = 1, n
        Q(i) = B (i)/ A (i, i)
    enddo 

    X = 0
    Xk = Q
        
    do while (sqrt(sum((Xk - X)**2)) >= eps)
        X = Xk
        do i = 1, n
            Xk(i) = dot_product(P(i, 1:i-1), Xk(1:i-1)) + dot_product(P(i, i+1:n), X(i+1:n)) + Q (i)
        enddo 
    enddo

    X = Xk
        
end function zeidel 


function relaxation (A, B) result (X)
    implicit none
    real :: A (1:, 1:), B (1:)
    real :: X (1:size(B))
    real :: P(1:size(B), 1:size(B))
    real :: Q(1:size(B)), Qk(1:size(B))
    integer :: n, i, j

    do i = 1, n 
        if(abs(A(i, i)) < sum(abs(A(i,1:i-1))) + sum(abs(A(i,i+1:n)))) then       ! проверка на диагональное преобладаение
            write(*, "('Матрица не обладает диагональным преобладанием')")
        endif
    enddo

    n = size (B)

    do i=1,n
        do j=1,n
            P(i,j) = - A(i,j)/A(i,i)
        enddo
    enddo

    do i = 1, n
        Q(i) = B (i)/ A (i, i)
    enddo 
        
    X = 0
    Qk = Q 

    do while (maxval(abs(Q)) >= eps)
        Q = Qk
        j = maxloc(abs(Q), dim = 1)
        X(j) = X(j) + Q(j)        
        
        do i = 1, n 
            Qk(i) = Q(i) + P(i, j)*Q(j)
        enddo 
    enddo  
end function relaxation


end module iterative_methods_sle
