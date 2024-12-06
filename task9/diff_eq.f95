module diff_eq
	use :: settings
	use :: gauss_integr
	use :: nonlinear_system
    implicit none

    contains

function runge_kutta(f, t1, x0) result(res)

	real, intent(in) :: t1(:), x0(:)
    real :: h, res(size(x0), size(t1))
    real :: k1(size(x0)), k2(size(x0)), k3(size(x0)), k4(size(x0))
    integer :: n, d, i

    interface
            pure function f(t, x) result(x_dot)
                real, intent(in) :: t, x(:)
                real :: x_dot(size(x))
            end function
    end interface
	n = size(t1)
	d = size(x0)
	h = t1(2)-t1(1)
	res(:,1) = x0
	do i=2,n
		k1 = h*  f(t1(i-1), res(:,i-1))
		k2 = h * f(t1(i-1)+h/2, res(:,i-1)+k1/2)
        k3 = h * f(t1(i-1)+h/2, res(:,i-1)+k2/2)
        k4 = h * f(t1(i-1), res(:,i-1)+k3)
        res(:,i) = res(:,i-1) + (k1 + 2*k2 + 2*k3 + k4) / 6
    end do
end function runge_kutta


function adams_extrapolation(f, t, x0, m) result(x)
	 interface
            pure function f(t, x) result(x_dot)
                real, intent(in) :: t, x(:)
                real :: x_dot(size(x))
            end function
    end interface
        real, intent(in) :: t(:), x0(:)
        integer, intent(in) :: m 
        real :: h, x(size(x0), size(t)), y(size(x0), size(t)), a(m)
        integer :: n, d, i

        n = size(t) !количество шагов
        h = t(2)-t(1) 
        a = make_A(m)
        x(:,:m) = runge_kutta(f, t(:m), x0) ! запуск методом Рунге-Кутты
        do i=1,m-1
            y(:,i) = f(t(i), x(:,i))
        end do
        do i=m+1,n
            y(:,i-1) = f(t(i-1), x(:,i-1))
            x(:,i) = x(:,i-1) + h * matmul(y(:,i-1:i-m:-1), a)
        end do
	end function adams_extrapolation
function make_A(n) result(a)
			integer, intent(in) :: n
			integer :: j
			real :: a(0:n-1)
			real :: coeff
			coeff = 1
			do j=2,n-1
				coeff = coeff * real(j)
			enddo
			coeff = 1.0/ coeff
			a(0) = integral_adamsA(0,n) * coeff

			do j = 1, n-2
				coeff = -coeff * real(n-j)/ real(j)
				a(j) = coeff*integral_adamsA(j,n)
			enddo
			coeff = -coeff /real(n-1)
			a(n-1) = coeff* integral_adamsA(n-1,n)
			end function make_A    

function integral_adamsA(j, n) result(res) 
		integer, intent(in) :: j, n
        real :: res
        res = gauss_int(fadamA, 0.0, 1.0, 8)
        contains
        pure function fadamA(z) result(res)
            real, intent(in) :: z
            real :: res
            integer :: k
            res = 1.0  ! Начальное значение для произведения
				do k = 0, n-1
					res = res * (z + real(k))
				end do
			res = res / (z + real(j))
        end function fadamA
end function integral_adamsA


function adams_interpolation(f, t, x0, m) result(x)
	 interface
            pure function f(t, x) result(x_dot)
                real, intent(in) :: t, x(:)
                real :: x_dot(size(x))
            end function
    end interface
        real, intent(in) :: t(:), x0(:)
        integer, intent(in) :: m 
        real :: h, x(size(x0), size(t)), y(size(x0), size(t)), b(m)
        integer :: n, d, i

        n = size(t) !количество шагов
        d = size(x0)
        h = t(2)-t(1) 
        b = make_B(m)
        x(:,:m) = runge_kutta(f, t(:m), x0) ! запуск методом Рунге-Кутты
        do i=1,m-1
            y(:,i) = f(t(i), x(:,i))
        end do
         do i=m,n 
            x(:,i) = newton(nonlinear, x(:,i-1), 100)
            y(:,i) = f(t(i), x(:,i))
        end do
	contains
		pure function nonlinear(XX) result(res)
            real, intent(in) :: XX(:)
            real :: y_tmp(d,m), res(d)
            y_tmp = y(:,i:i-m+1:-1)
            y_tmp(:,1) = f(t(i), XX)
            res = XX - x(:,i-1) - h * matmul(y_tmp, b)
        end function
        
       
	end function adams_interpolation
function make_B(n) result(b)
			integer, intent(in) :: n
			integer :: j
			real :: b(0:n-1)
			real :: coeff
			coeff = 1
			do j=1,n-1
				coeff = coeff * real(j)
			enddo
			coeff = 1.0/ coeff
			b(0) = integral_adamsB(-1,n) * coeff
			do j = 0, n-3
				coeff = -coeff * real(n-1-j)/ real(j+1)
				b(j+1) = coeff*integral_adamsB(j,n)
			enddo
			coeff = -coeff *real(n-3)/ real(n-1)
			b(n-1) = coeff* integral_adamsB(n-1,n)
			end function make_B    

function integral_adamsB(j, n) result(res) 
		integer, intent(in) :: j, n
        real :: res
        res = gauss_int(fadamB, 0.0, 1.0, 9)
        contains
        pure function fadamB(z) result(res)
            real, intent(in) :: z
            real :: res
            integer :: k
            res = 1.0  ! Начальное значение для произведения
				do k = -1, n-2
					res = res * (z + real(k))
				end do
			res = res / (z + real(j))
        end function fadamB
end function integral_adamsB   
                

end module diff_eq
    
