This repository contains some scripts on signal processing that could be useful for beginners who are learning those topics and want to see how the code works

numerical_ode.py is the solution of a ordinary differential equation using numerical methods, it is essentially a simulation of how heat flows in a metal plate, using the Heat Equation. The first function does the LU decomposition of a matrix,
the second funtion solves a linear system of a upper triangular matrix and main function solves the ODE https://en.wikipedia.org/wiki/LU_decomposition https://en.wikipedia.org/wiki/Heat_equation 

adaptative_filter.py has the implematation of the LMS and NLMS algorithms, that are the basis for projetcs in active noise cancelling https://en.wikipedia.org/wiki/Adaptive_filter

etimacaoespectral.py has a algorithm that gets spectral information from a signal, it uses the FFT as a basis but it does that in a form that works with stochastic data https://en.wikipedia.org/wiki/Fast_Fourier_transform#

kalman_filter.py is an implematation of a Kalman filter, the input "x" is a random walk evolving in time, the measured output "y" is the input "x" passed by a one parameter finite impulse response filter with added measurement noise https://en.wikipedia.org/wiki/Kalman_filter#

The code only uses numpy and matplotlib as libraries.
