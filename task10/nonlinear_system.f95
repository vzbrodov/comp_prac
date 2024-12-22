module nonlinear_system
	use :: my_prec
    use :: linear_system
    implicit none
    real(mp), parameter :: epsilon = 1e-7_mp
    contains

function newton (f, x0, iter) result (x)
    interface
        function f (y) result (fi)
			use :: my_prec
            real(mp), intent(in),dimension (1:) :: y
            real(mp), dimension (1:size(y)) :: fi
        end function f
    end interface

    real(mp), dimension (1:) :: x0
    integer :: i, iter
    real(mp), dimension (1:size(x0)) :: x, eps
    real(mp), dimension(1:size(x0), 1:size(x0)) :: H

    x = x0

    do i = 1, iter
        H = matrix_jacobi(f, x) ! Нахождение матрицы частных производных
        x = choice (H, matmul(H, x) - f(x)) ! Схема с выбором ведущего элемента
        eps = f(x) - f(x0)
        if (dot_product(eps, eps) < epsilon) return 
    enddo

end function newton

function matrix_jacobi(f, x) result (J)
    interface
        function f (y) result (fi)
			use :: my_prec
            real(mp), intent(in),dimension (1:) :: y
            real(mp), dimension (1:size(y)) :: fi
        end function f
    end interface
    real(mp), dimension (1:) :: x
    real(mp), dimension (1:size(x), 1:size(x)) :: J
    real(mp), dimension (1:size(x)) :: h
    real(mp), parameter :: sqreps = 1e-5_mp
    integer :: i
	
    do  i=1, size (x)
        h = 0.0_mp
        h(i) = sqreps
        J(:, i) = (f(x+h) - f(x-h))/(2.0_mp*h(i))
    end do
end function matrix_jacobi

end module nonlinear_system
    
