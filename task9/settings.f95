module settings
    implicit none

    integer, parameter :: end = 10
    real, dimension(2), parameter :: x0 = [1.0, -2.0]
	real, parameter :: h = 0.1
	integer, parameter :: k = 4 ! порядок для методов Адамса
contains 
	pure function f(t, x) result(x_dot)
		real, intent(in) :: t, x(:)
		real :: x_dot(size(x))
        x_dot(1) =  x(1) - x(2)
        x_dot(2) =  x(1) - x(2) - 2.0
    end function
    
end module settings
