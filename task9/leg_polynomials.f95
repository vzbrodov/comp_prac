module leg_polynomials
    implicit none
    real, parameter :: eps = 1e-7
    contains

function horner(a, x0) result(b) !деление на x-x0
    implicit none
        real, intent(in) :: a(0:), x0
        real :: b(0:size(a)-1)
        integer :: n, i
        n = size(a) - 1

        b(0) = a(0)
        do i = 1,n-1
            b(i) = b(i-1) * x0 + a(i)
        end do
end function

function make_legendre (n) result (p)
    implicit none
    integer,intent(in) :: n
    real, dimension (0:n) :: p
    integer :: i
    real, dimension (0:n) :: p1, p2 ! p1=p_(n-1), p2=p_(n-2)

    p1 = 0.0
    p2 = 0.0

    select case (n)
        case (0)
            p = 1.0
        case (1)
            p(0) = 1.0
            p(1) = 0.0
        case default
            p2(n) = 1.0 ! p_0(x) = 1
            p1(n-1) = 1.0 ! p_1(x) = x

            do i = 2,n
                p(:n-1) = p1(1:)
                p(n) = 0.0
                p =((2.0*real(i)- 1.0) * p - (real(i)-1.0) * p2) / real(i)
                p2 = p1
                p1 = p
            end do
    end select
end function make_legendre

function leg_roots(p) result(x)
    real, intent(in) :: p(0:)
    real :: x(1:size(p)-1), a1(0:size(p)-1)
    integer :: n, k
    n = size(p) - 1

    select case (n)
        case (1)
            x = - p(1) / p(0)
        case (2)
            x(1) = sqrt(-p(2)/p(0))
            x(2) = -sqrt(-p(2)/p(0))
        case (3:)
			a1 = p
			if ( mod(n,2) == 0 ) then
				do k = n,4,-2
					x(k) = sqrt(max_root(a1(0:k:2)))
					x(k-1)= -x(k)	
					a1(0:k-1) = horner(a1(0:k), x(k))
					a1(0:k-2) = horner(a1(0:k-1), x(k-1))
				end do
				x(1) = -sqrt(-a1(2)/a1(0))
				x(2) = sqrt(-a1(2)/a1(0))
			else 
				! a1 = horner(a1, real(0,16))
				x(1) = 0
				do k = n,3,-2
					if (k == 3) exit
					x(k) = sqrt(max_root(a1(0:k:2)))
					x(k-1)= -x(k)
					a1(0:k-1) = horner(a1(0:k), x(k))
					a1(0:k-2) = horner(a1(0:k-1), x(k-1))
				end do
				x(2) = -sqrt(-a1(2)/a1(0))
				x(3) = sqrt(-a1(2)/a1(0))
			endif
    end select

end function leg_roots

function max_root(a) result(x1)
    real, intent(in) :: a(0:)
    integer, parameter :: MAX_ITER = 1e6
    real :: x1, y(1:size(a)-1), b(1:size(a)-1)
    integer :: n, i

    n = size(a) - 1     
	b = a(1:) / a(0) 
	y = 1
	do i=1, MAX_ITER
		x1 = - dot_product(y, b)
		y(2:n) = y(1:n-1) 
		y(1) = x1
		x1 = - dot_product(y, b)
		y(2:n) = y(1:n-1) 
		y(1) = x1
		do while (abs(y(1)) < 0.03125)
                y = y * 1024
        end do
 
        x1 = y(1)/y(2)
		if (abs(x1 - y(2)/y(3)) < eps) exit
	enddo
end function max_root

end module leg_polynomials
