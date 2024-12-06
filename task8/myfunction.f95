module myfunction
    implicit none
    contains

function f(x) result (y)
implicit none
    real, intent(in) :: x
    real :: y
    y = x*x
end function f 

end module myfunction
