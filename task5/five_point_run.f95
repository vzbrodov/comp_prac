module five_point_run
implicit none

contains


function fivepoint_run_method(a0, b0) result (X)
    implicit none
    real, intent(in) :: a0(1:, 1:), b0(1:)
    real :: X(1:size(b0))
    real, dimension (-1:size(b0)) :: a,b,c !величины с индексами i<1 равны нулю, смещаю с 1 на -1
    real, dimension (-1:size(b0)) :: beta, alpha, p, q, r
    integer :: i, j, n
    
    n = size (b0)
    a = 0
    b = 0
    c = 0
    a(1:) = a0(3,:)
    b(1:) = a0(4,:)    
    c(1:) = a0(5,:) 
    p = 0
    q = 0
    r = 0
        do i = 1,n
            beta(i) = b(i-1) - p(i-2) * c(i-2)
            alpha(i) = a(i) - p(i-1) * beta(i) - q(i-2) * c(i-2)
            p(i) = (b(i) - q(i-1) * beta(i)) / alpha(i)
            q(i) = c(i) / alpha(i)
            r(i) = (b0(i) - r(i-1) * beta(i) - r(i-2) * c(i-2)) / alpha(i)
        end do
        x(n) = r(n)
        x(n-1) = r(n-1) - p(n-1) * x(n)
        do i = n-2,1,-1
            x(i) = r(i) - p(i) * x(i+1) - q(i) * x(i+2)
        end do

end function fivepoint_run_method

end module five_point_run