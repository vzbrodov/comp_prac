program main
    use :: fast_fourier_transform
implicit none
    integer, parameter :: IN = 1, OUT = 2, ABSf = 3
    integer :: n, i
    complex, allocatable, dimension (:) :: array, result
    real :: re, im
    
    open (IN, file = 'in/data.dat')        
    read (IN, '(2x, i10)') n 
        if (2**floor(log(real(n))/log(2.0)) .ne. n) then
            n = 2**(floor(log(real(n))/log(2.0)) + 1)
        endif
    
        allocate (array(0:n-1))
        array = 0

        do i=0,n-1
            read(IN, *) re, im
            array(i) = complex(re, im)
        enddo
    close (IN)  

    result = recFFT(array, -1) /sqrt(real(n))
    !result = recFFT(result, 1)/sqrt(real(n))
    !result = recFFT(recFFT(array, -1) /sqrt(real(n)),1) /sqrt(real(n))
    open (OUT, file = "out/result.dat")
    open (ABSf, file = "out/abs.dat")
        write(OUT, "('# ', i10)") size(result)
        do i = 0, size(result)-1
            write(ABSf, *) abs(result(i))
            write(OUT, *) real(result(i)), aimag(result(i))
        enddo
    close(OUT)
    close(ABSf)
deallocate(array, result)

end program main
