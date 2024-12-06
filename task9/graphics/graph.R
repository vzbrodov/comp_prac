t = seq(0,10,0.1)

plot(rk$X1,rk$X2, type = 'p', col='blue', xlab='t',ylab='X(i)')
title(main = "Runge_Kutta method")
lines(rk$X1,rk$X3,col="green",type = 'p')
lines(t,t^2+3*t+1, col='red')
lines(t,t^2+t-2, col='black')
legend("topleft", legend = c("x(1)_found", "x(2)_found","x(1)_analytic","x(2)_analytic"),
       lwd = 4, col =c("blue", "green","red","black"))


plot(ai$X1,ai$X2, type = 'p', col='blue', xlab='t',ylab='X(i)')
title(main = "Adams interpolation method ")
lines(ai$X1,ai$X3,col="green",type = 'p')
lines(t,t^2+3*t+1, col='red')
lines(t,t^2+t-2, col='black')
legend("topleft", legend = c("x(1)_found", "x(2)_found","x(1)_analytic","x(2)_analytic"),
       lwd = 4, col =c("blue", "green","red","black"))

plot(ae$X1,ae$X2, type = 'p', col='blue', xlab='t',ylab='X(i)')
title(main = "Adams extrapolation method")
lines(ae$X1,ae$X3,col="green",type = 'p')
lines(t,t^2+3*t+1, col='red')
lines(t,t^2+t-2, col='black')
legend("topleft", legend = c("x(1)_found", "x(2)_found","x(1)_analytic","x(2)_analytic"),
       lwd = 4, col =c("blue", "green","red","black"))

