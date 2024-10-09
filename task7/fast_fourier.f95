module fast_fourier_transform

    implicit none
        real, parameter :: pi = 4.0*atan(1.0)
    contains

recursive function recFFT(array0, sign) result (res1)
    complex, intent(in)  :: array0(0:)
    integer, intent(in)  :: sign
    complex :: res1(0:size(array0)-1), w(0:size(array0)/2 -1 )
    integer :: n, i

    n = size(array0)
  
    if (n .eq. 2) then
        res1(0) = array0(0) + array0(1)
        res1(1) = array0(0) - array0(1)
    else
        do i=0, n/2 - 1
            w(i) = exp(sign * cmplx(0,2) * pi * i / n)
        enddo
        res1(0::2) = recFFT(array0(:n/2-1) + array0(n/2:), sign)
        res1(1::2) = recFFT(w * (array0(:n/2-1) - array0(n/2:)), sign)
    endif

end function recFFT

end module fast_fourier_transform
    