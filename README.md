# Dynamics
This repository contains code which i created, used and updated before I get my Master's degree.

Studienarbeit came from a short essay in which I presented features of micro-slip model, including analytical one-bar model, discrete one-bar model and the brush model. Analytical one-bar model is derivated with the instructions from this book:

Title	Modelling Microslip Friction Damping and Its Influence on Turbine Blade Vibrations
Volume 519 of Linköping studies in science and technology: Dissertations, ISSN 0345-7524
Author	Gabor Csaba
Contributor	Universitetet i Linköping. Department of Mechanical Engineering. Division of Machine Design
Publisher	Division of Machine Design, Department of Mechanical Engineering, Linköping University, 1998
ISBN	9172191686, 9789172191686
Length	216 pages

Many thanks to Professor Csaba.

The other two models use numerical methods and is much simpler to understand.


DERIVEST SUIT is a set of numerical jacobian calculation tool. It has its own readme file and it works just fine. Slow however, but numerical.


Masterarbeit contains code used in my master thesis. The model is two masses on the ground connected with each other and fixed point through springs. The mass 1 has friction damping with ground and is bounded with fixed boundaries. The second mass is bounded with the first one with spring and bears the load, it has no frictional contact with ground. Both masses has viscous damping. The load is pre-defined. The original model is one-dimensional model which moves only in x-direction. The two-dimensional models(2mal 1d or real 2d) move in the xy-plane but their contact interfaces are different.

The contact interface of 2mal1d model is a simple extension of which from 1d model: the first mass moves on the ground and bears frictional forces in both directions. These forces are determined individually through different coefficient and displacement in x-and y-directions separately. They have no influence on each other. 

The contact interface of 2mal1d model is a simple extension of which from 1d model: the first mass moves on the ground and bears frictional forces only in one direction: the synthetic movement vector on the xy-plane. And there is only one friction coefficient in this case. 

Both models are solved in frequency domain which demands complex mass-and stiffniss matrices and dynamic load on the first mass represented in the form of fourier series. We consider however, in these programs only the first order from fourier series is taken.

The contact forces need to be calculated in time domain so an Alternating Frequency-Time method(AFT)is used. Only the coefficients of fourier transformation of nonlinear contact forces are extracted and used in the Balance equation，which is defined in the frequency domain and solved through Newton iteration. At each iteration the coefficients from Balance equation would be passed to the subroutine which transmitts it to time domain with IFFT and updates it through slip-stick transition in time domain and provides an updated complex Jacobian matrix in the end for Newton's method in frequency domain.

