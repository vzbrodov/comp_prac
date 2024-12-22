module settings
	use :: my_prec
    implicit none

    real, parameter :: t_end = 50.0_mp
    real(mp), dimension(2), parameter :: x0 = [1.5_mp, 1.0_mp]
	real(mp), parameter :: h = 0.01_mp
	real(mp), parameter :: mu = 10.0_mp
	integer, parameter :: k = 6 ! порядок для методов Адамса
contains 
	pure function f( x) result(x_dot)
		real(mp), intent(in) ::  x(:) 
		real(mp) :: x_dot(size(x))
        x_dot(1) = mu*(x(1) - x(1)*x(1)*x(1)/3.0_mp - x(2))
		x_dot(2) = x(1) / mu
    end function
    
end module settings
