module tridiagonal
implicit none
contains

function tridiag (A, B) result (C)
implicit none

        real, dimension (1:,1:) :: A, B
        real, dimension (size(A(:,1)),1:5) :: C

        integer n, i, j, k 

        n = size (A(:,1))

                !Первые две и две последние строки в произведении отличаются от остальных строк количеством ненулевых элементов.
                !можно посчитать отдельно три блока: первые две, последние две и строки со второй до n-3 
                C(1,1) = A(1,1)*B(1,1) + A(1,2)*B(2,1)      
                C(1,2) = A(1,1)*B(1,2) + A(1,2)*B(2,2)
                C(1,3) = A(1,2)*B(2,3)
                C(1,4) = 0                                      
                C(1,5) = 0

                C(2,1) = A(2,1)*B(1,1) + A(2,2)*B(2,1)        
                C(2,2) = A(2,1)*B(1,2) + A(2,2)*B(2,2) + A(2,3)*B(3,1)
                C(2,3) = A(2,2)*B(2,3) + A(2,3)*B(3,2) 
                C(2,4) = A(2,3)*B(3,3)
                C(2,5) = 0

        do i = 2, n - 3   
                                       !строчки с 3 до n-2 (из за того, что матрицы А и B вида (n x 3), делаем "почти" умножение матриц, но только 
                k = i + 1                !со смещенными коэффциентами на 1 и на 2
                j = i + 2

                C(k,1) = A(k,1)*B(i,1)
                C(k,2) = A(k,1)*B(i,2) + A(k,2)*B(k,1)
                C(k,3) = A(k,1)*B(i,3) + A(k,2)*B(k,2) + A(k,3)*B(j,1) 
                C(k,4) = A(k,2)*B(k,3) + A(k,3)*B(j,2)
                C(k,5) = A(k,3)*B(j,3)
        enddo

                i = n - 2
                k = n - 1 
                j = n                 

                C(k,1) = 0                    
                C(k,2) = A(k,1)*B(i,1) 
                C(k,3) = A(k,1)*B(i,2) + A(k,2)*B(k,1)
                C(k,4) = A(k,1)*B(i,3) + A(k,2)*B(k,2) + A(k,3)*B(j,2) 
                C(k,5) = A(k,2)*B(k,3) + A(k,3)*B(j,3)

                C(j,1) = 0                      
                C(j,2) = 0
                C(j,3) = A(j,2)*B(k,1)
                C(j,4) = A(j,2)*B(k,2) + A(j,3)*B(j,2)
                C(j,5) = A(j,2)*B(k,3) + A(j,3)*B(j,3)

end function tridiag
end module tridiagonal
