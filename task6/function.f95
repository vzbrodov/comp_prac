module myfunction
	implicit none
	contains

function f(x) result (fi)
implicit none
	real, dimension (1:) :: x
    real, dimension (1:size(x)) :: fi

 	fi(1)=x(1) - x(2)**2 
    fi(2) = cos(x(1)) - x(2)
end function f 

end module myfunction