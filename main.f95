program one
        use omp_lib
        implicit none

        integer n, i, j, k
        real(16), dimension (:,:), allocatable :: A, B, C
        integer, parameter :: data1 = 1, data2 = 2
                open (data1, file = './data1.dat')
                open (data2, file = './data2.dat')

                        read (data1, '(2x, i10)') n
                        read (data2, '()')

                        allocate (A(n, n), B(n, n), C(n, n))

                        read(data1, *) A    
                        read(data2, *) B     

                        A = transpose (A)
                        B = transpose (B)


                close(data1)
                close(data2)

        !$omp parallel do 

                do i = 1, n
                        do j = 1, n
                                C(i, j) = 0
                                do k = 1, n
                                        C(i, j) = C(i, j) + A(i, k) * B(k, j)
                                enddo
                        enddo
                enddo

        !$omp end parallel do

                open (3, file = 'result.dat')
                        write(3, '("# "i10)') n

                        do i = 1, n
                                write(3, '(10000(f16.5),3x)') (C(i, j), j = 1, n)
                        end do
                close (3)
        deallocate (A, B, C)


end program one

