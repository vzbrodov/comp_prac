module settings
	use :: my_prec
    implicit none

    real, parameter :: t_end = 1e-6_mp
    real(mp), dimension(3), parameter :: x0 = [0.1_mp, 0.1_mp, 0.1_mp]
	real(mp), parameter :: h = 1.0e-9_mp
	integer, parameter :: k = 10 ! порядок для методов Адамса
contains 
	pure function f( x) result(x_dot)
		real(mp), intent(in) ::  x(:) 
		real(mp) :: x_dot(size(x))
        x_dot(1) = -0.05_mp*x(1) + 1e4_mp*x(2)*x(3)
		x_dot(2) = 0.05_mp*x(1) - 1e4_mp*x(2)*x(3) - 1e7_mp*x(2)
		x_dot(3) = 1e7_mp*x(2)
    end function
    
end module settings
