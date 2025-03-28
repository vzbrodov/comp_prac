module linear_system
    implicit none
    real, parameter :: eps = 1e-3
    contains

function gauss (A,B) result (X)
    implicit none
    real :: A(1:, 1:), B(1:)
    real :: X(1:size(B))
    integer :: i, j, n

    n = size (B)

    do i = 1, n !приводим матрицу к верхней треугольной
        if (abs(A(i,i)) < eps) then
            write(*, *)'Деление на близкое к нулю число' 
        endif

        B(i) = B(i)/A(i,i)    
        A(i,:) = A(i,:)/A(i,i)

            do j = i+1, n
                B(j) = B(j) - B(i)*A(j,i)
                A(j,:) = A(j,:) - A(i,:)*A(j,i)
            enddo 
    enddo

    X(n) = B(n) !после этого шага матрица верхняя треугольная
    !обратный ход метода гаусса 
    do i = n-1, 1, -1
        X(i) = B(i)
        do j = i+1, n
            X(i) = X(i) - A(i,j) * X(j)
        end do
    enddo
end function gauss


function jordan (A, B) result (X)
    implicit none
    real :: A(1:,1:), B(1:)
    real :: X(1:size(B))
    integer :: i, j, n

    n = size (B)

    do i = 1, n
        if (abs(A(i,i)) < eps) then
            write(*,*)'Деление на близкое к нулю число' 
        endif

        B(i) = B(i) / A(i,i)       
        A(i, :) = A(i,:) / A(i,i)
        do j = 1, n
            if (j == i) cycle                
            B(j) = B(j) - B(i)*A(j,i)
            A(j,:) = A(j,:) - A(i,:)*A(j,i)
        enddo 
    enddo
    X = B

end function jordan 

function choice (A,B) result (X)
    real :: A(1:, 1:), B(1:)
    real :: X(1:size(B))
    integer :: n, i, j
    integer :: idx(2), transp(2, size(B))

    n = size (B)

    do i = 1, n
        idx = maxloc(abs(A(i:n,i:n))) + i - 1 ! поиск индексов максимального элемента
            call swap_vect(A(i,:),A(idx(1),:)) ! swap_vect(row1,row2) - меняет строки местами 
            call swap_sc(B(i),B(idx(1)))        ! swap_sc(b,a) - меняет скалярные значения местами 
            call swap_vect(A(:,i), A(:,idx(2)))
        transp(:,i) = idx   !записываем в каком порядке перестановка строк

        if (abs(A(i,i)) < eps) then
            write(*,*)'Деление на близкое к нулю число' 
        endif

        B(i) = B(i) / A(i,i)       
        A(i,:) = A(i,:) / A(i,i)
        do j = 1, n
            if (j == i) cycle                
            B(j) = B(j) - B(i)*A(j,i)
            A(j,:) = A(j,:) - A(i,:)*A(j,i)
        enddo 
    enddo

    do i = n, 1, -1 ! возвращаем правильный порядок строк
        idx = transp(:,i)
            call swap_vect(A(i,:),A(idx(1),:))
            call swap_sc(B(i),B(idx(1))) 
            call swap_vect(A(:,i), A(:,idx(2)))
        enddo

    do i = 1, n
        X(maxloc(A(i,:))) = B(i)
    enddo

end function choice


subroutine swap_vect(row1,row2)
    real :: row1(:), row2(:)
    real :: tmp
    integer :: i
    do i = 1, size(row1)
        tmp = row1(i)
        row1(i) = row2(i)
        row2(i) = tmp
    end do
end subroutine swap_vect

subroutine swap_sc(a,b)
        real :: a, b
        real :: tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap_sc
end module linear_system
    