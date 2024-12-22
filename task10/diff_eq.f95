module diff_eq
	use :: my_prec
	use :: settings
	use :: gauss_integr
	use :: nonlinear_system
	use :: linear_system
    implicit none

    contains

function runge_kutta(f, t1, x0) result(res)
	real(mp), intent(in) :: t1(:), x0(:)
    real(mp) :: h, res(size(x0), size(t1))
    real(mp) :: k1(size(x0)), k2(size(x0)), k3(size(x0)), k4(size(x0))
    integer :: n, d, i

    interface
		pure function f(x) result(x_dot)
			use :: my_prec
            real(mp), intent(in) :: x(:)
            real(mp) :: x_dot(size(x))
        end function
    end interface
	n = size(t1)
	d = size(x0)
	h = t1(2)-t1(1)
	res(:,1) = x0

    do i = 1, n-1
		k1 = h*f(res(:,i ))
        k2 = h*f(res(:,i) + k1/2.0_mp)
        k3 = h*f(res(:, i) + k2/2.0_mp)
        k4 = h*f(res(:,i) + k3)
        res(:, i+1) = res( :, i) + (k1 + 2.0_mp*k2 + 2.0_mp*k3 + k4)/6.0_mp
    enddo
end function runge_kutta

subroutine adams_extrapolation(f, t, x0, a, m, x_pred)
  interface
    pure function f(x) result(x_dot)
      use :: my_prec
      real(mp), intent(in) :: x(:)
      real(mp) :: x_dot(size(x))
    end function
  end interface
  real(mp), intent(in) :: t(:), x0(:), a(:)
  integer, intent(in) :: m
  real(mp), intent(out) :: x_pred(size(x0))
  real(mp) :: h, x(size(x0), size(t)), y(size(x0), size(t))
  integer :: n, d, i

  n = size(t) ! количество шагов
  h = t(2) - t(1)
  x(:,:m) = runge_kutta(f, t(:m), x0) ! запуск методом Рунге-Кутты
  do i = 1, m - 1
    y(:,i) = f(x(:,i))
  end do
  do i = m + 1, n
    y(:,i-1) = f(x(:,i-1))
    x_pred = x(:,i-1) + h * matmul(y(:,i-1:i-m:-1), a)
  end do
end subroutine adams_extrapolation

function adams_interpolation_from9(f, t, x0, m) result(x)
	 interface
    pure function f(x) result(x_dot)
      use :: my_prec
      real(mp), intent(in) :: x(:)
      real(mp) :: x_dot(size(x))
    end function
  end interface
        real(mp), intent(in) :: t(:), x0(:)
        integer, intent(in) :: m 
        real(mp) :: h, x(size(x0), size(t)), y(size(x0), size(t)), b(m)
        integer :: n, d, i

        n = size(t) !количество шагов
        d = size(x0)
        h = t(2)-t(1) 
        b = make_B(m)
        x(:,:m) = runge_kutta(f, t(:m), x0) ! запуск методом Рунге-Кутты
        do i=1,m-1
            y(:,i) = f(x(:,i))
        end do
         do i=m,n 
            x(:,i) = newton(nonlinear, x(:,i-1), 100)
            y(:,i) = f( x(:,i))
        end do
	contains
		pure function nonlinear(XX) result(res)
            real(mp), intent(in) :: XX(:)
            real(mp) :: y_tmp(d,m), res(d)
            y_tmp = y(:,i:i-m+1:-1)
            y_tmp(:,1) = f(XX)
            res = XX - x(:,i-1) - h * matmul(y_tmp, b)
        end function
        
       
	end function adams_interpolation_from9

subroutine adams_interpolation(f, t, x0, b, m, x_corr) !это модифицированный Адамс, для 10 задания, выводит только последний шаг
  interface
    pure function f(x) result(x_dot)
      use :: my_prec
      real(mp), intent(in) :: x(:)
      real(mp) :: x_dot(size(x))
    end function
  end interface
  real(mp), intent(in) :: t(:), x0(:), b(:)
  integer, intent(in) :: m
  real(mp), intent(out) :: x_corr(size(x0))
  real(mp) :: h, x(size(x0), size(t)), y(size(x0), size(t))
  integer :: n, d, i

  n = size(t) ! количество шагов
  d = size(x0)
  h = t(2) - t(1)
  x(:,:m) = runge_kutta(f, t(:m), x0) ! запуск методом Рунге-Кутты
  do i = 1, m - 1
    y(:,i) = f(x(:,i))
  end do
  do i = m, n
    x_corr = newton(nonlinear, x(:,i-1), 10)
    y(:,i) = f(x_corr)
  end do

contains
  pure function nonlinear(XX) result(res)
    use :: my_prec
    real(mp), intent(in) :: XX(:)
    real(mp) :: y_tmp(d,m), res(d)
    integer :: j


    do j = 1, m
      y_tmp(:,j) = y(:,i-j+1)
    end do
    y_tmp(:,1) = f(XX)
    res = XX - x(:,i-1) - h * matmul(y_tmp, b)
  end function nonlinear
end subroutine adams_interpolation

function make_A(n) result(a)
	integer, intent(in) :: n
	integer :: j
	real(mp) :: a(0:n-1)
	real(mp) :: coeff
	coeff = 1
	do j=2,n-1
		coeff = coeff * real(j, mp)
	enddo
	coeff = 1.0_mp/ coeff
	a(0) = integral_adamsA(0,n) * coeff
	do j = 1, n-2
		coeff = -coeff * real(n-j, mp)/ real(j, mp)
		a(j) = coeff*integral_adamsA(j,n)
	enddo
	coeff = -coeff /real(n-1,mp)
	a(n-1) = coeff* integral_adamsA(n-1,n)
end function make_A    

function integral_adamsA(j, n) result(res) 
	integer, intent(in) :: j, n
    real(mp) :: res
    res = gauss_int(fadamA, 0.0_mp, 1.0_mp, 10)
    contains
		pure function fadamA(z) result(res)
			use :: my_prec
            real(mp), intent(in) :: z
            real(mp) :: res
            integer :: k
            res = 1.0_mp  ! Начальное значение для произведения
			do k = 0, n-1
				res = res * (z + real(k, mp))
			end do
			res = res / (z + real(j, mp))
        end function fadamA
end function integral_adamsA

function make_B(n) result(b)
	integer, intent(in) :: n
	integer :: j
	real(mp) :: b(-1:n-2)
	real(mp) :: coeff
	coeff = 1.0_mp
	do j=1,n-1
		coeff = coeff * real(j,mp)
	enddo
	coeff = 1.0_mp/ coeff
	b(-1) = integral_adamsB(-1,n) * coeff
	do j = 0, n-3
		coeff = -coeff * real(n-1-j, mp)/ real(j+1, mp)
		b(j+1) = coeff*integral_adamsB(j,n)
	enddo
	coeff = -coeff / real(n-1, mp)
	b(n-2) = coeff* integral_adamsB(n-2,n)
end function make_B    

function integral_adamsB(j, n) result(res) 
	integer, intent(in) :: j, n
    real(mp) :: res
    res = gauss_int(fadamB, 0.0_mp, 1.0_mp, 12)
contains
	pure function fadamB(z) result(res)
		use :: my_prec
        real(mp), intent(in) :: z
        real(mp) :: res
        integer :: k
        res = 1.0_mp
			do k = -1, n-2
				res = res * (z + real(k, mp))
			end do
		res = res / (z + real(j, mp))
	end function fadamB
end function integral_adamsB   
                
function rosenbroke(f, t1, x0) result(res)
	real(mp), intent(in) :: t1(:), x0(:)
    real(mp) :: res(size(x0), size(t1)), h
    real(mp), dimension(size(x0),size(x0)) :: E, J, T
    integer :: n, d, i
	real(mp), parameter :: my_alpha = 1.077_mp
    real(mp), parameter :: my_beta = -0.372_mp
    real(mp), parameter :: my_gamma = -0.577_mp
    interface
		pure function f(x) result(x_dot)
			use :: my_prec	
			real(mp), intent(in) ::  x(:)
            real(mp) :: x_dot(size(x))
        end function
    end interface
    n = size(t1)
	d = size(x0)
	h = t1(2)-t1(1)
	J = matrix_jacobi(f, x0)
	E = 0.0_mp
	do i = 1, d
		E(i, i) = 1.0_mp
    enddo
	T = E - my_alpha*h*J - my_beta*h*h*matmul(J, J)
	res(:,1) = x0
	do i = 1, n
		res(:, i+1) = choice(T, h*f(res(:, i) + my_gamma*h*f(res(:, i))) + matmul(T, res(:, i)))
    enddo
    
end function rosenbroke

function precor(f, t1, x0, m) result(x)
  interface
    pure function f(x) result(x_dot)
      use :: my_prec
      real(mp), intent(in) :: x(:)
      real(mp) :: x_dot(size(x))
    end function
  end interface
  real(mp), intent(in) :: x0(:), t1(:)
  real(mp) :: x(size(x0), size(t1)), a(m), b(m)
  integer, intent(in) :: m
  real(mp) :: h
  integer :: n, d, i, j
  real(mp), allocatable :: x_pred(:), x_corr(:)

  n = size(t1)
  d = size(x0)
  a = make_A(m)
  b = make_B(m)
  h = t1(2) - t1(1)
  x(:,:m) = runge_kutta(f, t1(:m), x0)
   allocate(x_pred(d), x_corr(d))
  do i = m, n


    call adams_extrapolation(f, t1(i-m+1:i), x(:,i-m+1), a, m, x_pred) 

    call adams_interpolation(f, t1(i-m+1:i), x(:,i-m+1), b, m, x_corr)
    
    x(:,i) = x_corr
    
  end do
  deallocate(x_pred, x_corr)
end function precor

end module diff_eq
    
