function [Fout] = calc_nlforce (qamp,kt,fn,mu)
% Calculates the first harmonic of the friction force
%kt,fn,mu,debug
%==========================================================================
% Inputs
%--------------------------------------------------------------------------
% q         = contact displacement
% kt        = tangential contact stiffness
% fnormal   = normal force
% mu        = coefficient of friction
%==========================================================================
% Outputs
%--------------------------------------------------------------------------
% f_nl     = harmonic contact force
%==========================================================================
% kt = 60e6;
% fn = 22;
% mu = 0.5;
% debug = 1;
% Exitation
qt = linspace(0,4*pi,2*2048+1); % don't change this, hard coded!
qt = qt(1:end-1);

qx = qamp(1)*cos(qt)+qamp(2)*sin(qt);
qy = qamp(3)*cos(qt)+qamp(4)*sin(qt);







N = 2*2048;

% Solving process

%initialization
x = zeros(length(qt),1);  % old Q
y = zeros(length(qt),1);  % old Q
ftx = zeros(length(qt),1);
fty = zeros(length(qt),1);

contactstatus = 0;
%========================================================================== 
% main iteration loop
%========================================================================== 
for it = 2: length(qt)
% ---------------Slip in all positive directions-------------------------------- 
    if kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn > 0 
        laenge = sqrt((qx(it)-x(it-1))^2+(qy(it)-y(it-1))^2);
        ftx(it) = -mu*fn*(qx(it)-x(it-1))/laenge;
        fty(it) = -mu*fn*(qy(it)-y(it-1))/laenge;        
        x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
        y(it) = qy(it) + fty(it)/kt;



% ---------------General stick---------------------------------------------
    else
        ftx(it) = -kt*(qx(it) - x(it-1));
        fty(it) = -kt*(qy(it) - y(it-1));
        x(it) = x(it-1);
        y(it) = y(it-1);      
    end

   
end

Fx = ffttho(ftx(length(ftx)/2+1:end),1,'nfft',length(ftx)/2); %only use second period for FFT
Fy = ffttho(fty(length(fty)/2+1:end),1,'nfft',length(fty)/2); %only use second period for FFT

Foutxreal  = 2*real(Fx(length(Fx)/2+2)); %first harmonic
Foutyreal  = 2*real(Fy(length(Fy)/2+2)); %first harmonic
Foutximg  = -2*imag(Fx(length(Fx)/2+2)); %first harmonic
Foutyimg  = -2*imag(Fy(length(Fy)/2+2)); %first harmonic

Fout = [Foutxreal;Foutximg;Foutyreal;Foutyimg]';
