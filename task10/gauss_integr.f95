module gauss_integr
    use :: leg_polynomials
    use :: linear_system
    use :: my_prec
    implicit none

contains

function gauss_int(f, a, b, n) result(res)
    real(mp), intent(in) :: a, b
    integer, intent(in) :: n
    integer :: i 
    character(10) :: filename
    logical :: file_exists
    real(mp) :: res, coeff(n), t(n), ft(n)
    interface
        function f(x) result(y)
			use :: my_prec
            real(mp), intent(in) :: x
            real(mp) :: y
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

	t = a + (t+1.0_mp) * (b-a) / 2.0_mp
    do i = 1, n
		ft(i) = f(t(i))
    enddo

    res = dot_product(coeff, ft) * (b-a) / 2.0_mp
end function

function gauss_quad_coeff(roots) result(A)
	real(mp), intent(in) :: roots(:)
	real(mp) :: t_matrix(1:size(roots),1:size(roots))
    real(mp) :: fr(size(roots)), A(size(roots))
    integer :: n, k
    n = size(roots)
        
	do k = 1, n
		t_matrix(k, :) = roots(:)**real(k-1)
        fr(k) = 2.0_mp*mod(k, 2)/ real(k)
	enddo
	A = gauss(t_matrix, fr)
end function

end module gauss_integr
