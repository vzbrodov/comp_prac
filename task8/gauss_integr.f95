module gauss_integr
    use :: leg_polynomials
    use :: linear_system
    implicit none

contains

function gauss_int(f, a, b, n) result(res)
    real, intent(in) :: a, b
    integer, intent(in) :: n
    integer :: i 
    character(10) :: filename
    logical :: file_exists
    real :: res, coeff(n), t(n), ft(n)
    interface
        function f(x) result(y)
            real, intent(in) :: x
            real :: y
        end function
    end interface
    write(filename, "('quad', i2.2, '.dat')") n


    INQUIRE(FILE=filename, EXIST=file_exists)

    if (file_exists) then
        open(1, file=filename)
            do i=1,n
                read(1,*) coeff(i), t(i)
            end do
        close(1)
    else
        t = leg_roots(make_legendre(n))
        coeff = gauss_quad_coeff(t)
        open(1, file=filename)
            do i=1,n
                write(1,*) coeff(i), t(i)
            end do
        close(1)
    endif

	t = a + (t+1.0) * (b-a) / 2.0
    do i = 1, n
		ft(i) = f(t(i))
    enddo

    res = dot_product(coeff, ft) * (b-a) / 2.0
end function

function gauss_quad_coeff(roots) result(A)
	real, intent(in) :: roots(:)
	real :: t_matrix(1:size(roots),1:size(roots))
    real :: fr(size(roots)), A(size(roots))
    integer :: n, k
    n = size(roots)
        
	do k = 1, n
		t_matrix(k, :) = roots(:)**real(k-1)
        fr(k) = 2.0*mod(k, 2)/ real(k)
	enddo
	A = gauss(t_matrix, fr)
end function

end module gauss_integr
