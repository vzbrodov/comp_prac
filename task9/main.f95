program main
    use :: settings
    use :: diff_eq
implicit none
	integer, parameter :: RungeKutta = 1
	integer, parameter :: extr_Adams = 2, inter_Adams = 3
	
	integer :: n, d, i

	real :: t(int(end/h)+1)
	real :: res1(size(x0),size(t))
	n = size(t)
    d = size(x0)
    
    do i=1,n
        t(i) = (i-1) * h
    end do
    
    res1 = runge_kutta(f, t, x0)
    open (RungeKutta, file = "rk.dat")
        do i = 1,n 
                write(RungeKutta, *) t(i), res1(:,i)
        enddo
	close(RungeKutta)
    res1 = adams_extrapolation(f, t, x0, k)
    open (extr_Adams, file = "ae.dat")
        do i = 1, n 
                write(extr_Adams, *) t(i), res1(:,i)
        enddo
	close(extr_Adams)
    
    res1 = adams_interpolation(f, t, x0, k)
    
    open (inter_Adams, file = "ai.dat")
        do i = 1,n 
                write(inter_Adams, *) t(i), res1(:,i)
        enddo
	close(inter_Adams)
end program main
