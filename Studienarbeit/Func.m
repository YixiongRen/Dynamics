function [ y ] = Func( x )
%FUNC Summary of this function goes here
%   Detailed explanation goes here

miu = 0.3;  %      empfindlich miu Inversely proportional to u
q0 = 100;           % Max_Ff
q2 = q0/2;           % Max_Ff    q Inversely proportional to u
E = 21000000000;        % beeinfluss die Verschiebung  E Inversely proportional to u
A = 0.0001;         % beeinfluss die Verschiebung  A Inversely proportional to u
L = 0.002;   %      empfindlich L proportional to u
F_amp = 50;         
Omega = 30;  
y = -4*miu*q2/3/L^2*x^3+2*miu*q2/L*x^2+miu*q0*x-F_amp;

end

