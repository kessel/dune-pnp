set term postscript enhanced eps color 

set out "one_wall.eps"
lambda=1./sqrt(8*pi)
kappa=1/lambda
F=.001
F=20
sigma=F/4/pi/4./pi/2
phi0=2*asinh(2*pi*kappa*sigma)
g(x)=2*log( (1+tanh(0.25*phi0)*exp(-kappa*x))/(1-tanh(0.25*phi0)*exp(-kappa*x)) )

plot "solution.dat.dat" u 1:(-$3), F/kappa*exp(-kappa*x), g(x)

!epstopdf one_wall.eps



