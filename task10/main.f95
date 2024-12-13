program main
    use :: settings
    use :: my_prec
    use :: diff_eq
implicit none
	integer, parameter :: Rosenbroke_f = 1
	integer, parameter :: Precor_f = 2
	
	integer :: n, d, i

	real(mp) :: t(int(t_end/h)+1)
	real(mp) :: res1(size(x0),size(t))
	n = size(t)
    d = size(x0)
    
    do i=1,n
        t(i) = (i-1) * h
    end do
    res1=0
    res1 = rosenbroke(f, t, x0)
    open (Rosenbroke_f, file = "rosen.dat")
        do i = 1,n 
                write(Rosenbroke_f, *) t(i), res1(:,i)
        enddo
	close(Rosenbroke_f)
	res1=0
    res1 = precor(f, t, x0, k)
    open (Precor_f, file = "precor.dat")
        do i = 1, n 
                write(Precor_f, *) t(i), res1(:,i)
        enddo
	close(Precor_f)
    
end program main
