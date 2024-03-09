program two
        use tridiagonal
        implicit none

        integer, parameter :: data1 = 1, data2 = 2
        integer n, i
        real, dimension (:,:), allocatable :: A, B, C 


        open (data1, file = './data1.dat')
        open (data2, file = './data2.dat')

        read (data1, '(2x, i10)') n
        read (data2, '()')

        allocate (A(n, 3), B(n, 3), C(n, 5))  !на входе две трехдиагональные матрицы, на выходе пятидиагональная.
        !можно трехдиагональные матрицы записывать как матрицы (n x 3), так как остальные будут нули, незачем их таскать и занимать память

        !так как в первой и последней строке одно элемента не существует(он равен ноль и мы его не пишем), прочитаем файл по отдельности
        !сначала 1 строку, со 2 по n-1, и последнюю
        read(data1, *) A(1, 1:2)
        read(data2, *) B(1, 1:2)
        
        do i=2,n-1
                read(data1, *) A(i, :)
                read(data2, *) B(i, :)
        enddo

        read(data1, *) A(n, 2:3)
        read(data2, *) B(n, 2:3)

        close (data1)
        close (data2)

                C = tridiag(A, B)
        open (3, file = 'result.dat')

                write (3, '("# "i10)') n
                write (3, '(3(f16.8, 1x))') C(1, 1), C(1, 2), C(1, 3)
                write (3, '(4(f16.8, 1x))') C(2, 1), C(2, 2), C(2, 3), C(2, 4)
                do i=3,n-2
                write (3, '(5(f16.8, 1x))') C(i,1:5)
                enddo 
                write (3, '(4(f16.8, 1x))') C(n-1, 2), C(n-1, 3), C(n-1, 4), C(n-1, 5)
                write (3, '(3(f16.8, 1x))') C(n, 3), C(n, 4), C(n, 5)

                close (3)


        deallocate (A, B, C)

end program two
