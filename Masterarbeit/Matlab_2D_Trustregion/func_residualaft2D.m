function [R J] = func_residualaft2D(qamp, S, fl, kt, fn, mu, debug)
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
dfxdrQoldx = 0;   % real part
dfxdiQoldx = 0;   % imaginary part
dfydrQoldy = 0;   % real part
dfydiQoldy = 0;   % imaginary part
dfydrQoldx = 0;   % real part
dfydiQoldx = 0;   % imaginary part
dfxdrQoldy = 0;   % real part
dfxdiQoldy = 0;   % imaginary part
dfxdrQx = zeros(length(qt),1);
dfxdiQx = zeros(length(qt),1);
dfydrQy = zeros(length(qt),1);
dfydiQy = zeros(length(qt),1);
dfxdrQy = zeros(length(qt),1);
dfxdiQy = zeros(length(qt),1);
dfydrQx = zeros(length(qt),1);
dfydiQx = zeros(length(qt),1);
statu = [0];
contactstatus = 0;
%========================================================================== 
% main iteration loop
%========================================================================== 
for it = 2: length(qt)
% ---------------General slip---------------------------------------------- 
    if kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn > 0
        laenge = sqrt((qx(it)-x(it-1))^2+(qy(it)-y(it-1))^2);
        ftx(it) = -mu*fn*(qx(it)-x(it-1))/laenge;
        fty(it) = -mu*fn*(qy(it)-y(it-1))/laenge;
        h = cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2));
        h_ = sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2));
        dfxdrQx(it) = -mu*fn*(1/laenge*(h-dfxdrQoldx/kt)-(qx(it)-x(it-1))/laenge^3*((qx(it)-x(it-1))*(h-dfxdrQoldx/kt)-(qy(it)-y(it-1))*dfydrQoldx/kt));
        dfxdiQx(it) = -mu*fn*(1/laenge*(h_-dfxdiQoldx/kt)-(qx(it)-x(it-1))/laenge^3*((qx(it)-x(it-1))*(h_-dfxdiQoldx/kt)-(qy(it)-y(it-1))*dfydiQoldx/kt));
        dfydrQy(it) = -mu*fn*(1/laenge*(h-dfydrQoldy/kt)-(qy(it)-y(it-1))/laenge^3*((qy(it)-y(it-1))*(h-dfydrQoldy/kt)-(qx(it)-x(it-1))*dfxdrQoldy/kt));
        dfydiQy(it) = -mu*fn*(1/laenge*(h_-dfydiQoldy/kt)-(qy(it)-y(it-1))/laenge^3*((qy(it)-y(it-1))*(h_-dfydiQoldy/kt)-(qx(it)-x(it-1))*dfxdiQoldy/kt));
        dfxdrQy(it) = -mu*fn*(-dfxdrQoldy/kt*(qy(it)-y(it-1))^2-(qx(it)-x(it-1))*(qy(it)-y(it-1))*(h-dfydrQoldy/kt))/laenge^3;
        dfxdiQy(it) = -mu*fn*(-dfxdiQoldy/kt*(qy(it)-y(it-1))^2-(qx(it)-x(it-1))*(qy(it)-y(it-1))*(h_-dfydiQoldy/kt))/laenge^3;
        dfydrQx(it) = -mu*fn*(-dfydrQoldx/kt*(qx(it)-x(it-1))^2-(qx(it)-x(it-1))*(qy(it)-y(it-1))*(h-dfxdrQoldx/kt))/laenge^3;
        dfydiQx(it) = -mu*fn*(-dfydiQoldx/kt*(qx(it)-x(it-1))^2-(qx(it)-x(it-1))*(qy(it)-y(it-1))*(h_-dfxdiQoldx/kt))/laenge^3;
        x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
        y(it) = qy(it) + fty(it)/kt;
        contactstatus = 11;
% ---------------General stick---------------------------------------------
    elseif kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn <= 0
        ftx(it) = -kt*(qx(it) - x(it-1));
        fty(it) = -kt*(qy(it) - y(it-1));
        dfxdrQx(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfxdrQoldx; % the factor 2*2*pi is hard coded due to 2 periods!
        dfxdiQx(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfxdiQoldx;
        dfydrQy(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfydrQoldy; % the factor 2*2*pi is hard coded due to 2 periods!
        dfydiQy(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfydiQoldy;
        dfxdrQy(it) = dfxdrQy(it-1);
        dfxdiQy(it) = dfxdiQy(it-1);
        dfydrQx(it) = dfydrQx(it-1);
        dfydiQx(it) = dfydiQx(it-1);
        x(it) = x(it-1);
        y(it) = y(it-1);
        contactstatus = 00;        
    end
    statu(end+1)=contactstatus;
    dfxdrQoldx = dfxdrQx(it);
    dfxdiQoldx = dfxdiQx(it);
    dfydrQoldy = dfydrQy(it);
    dfydiQoldy = dfydiQy(it);
    dfxdrQoldy = dfxdrQy(it);
    dfxdiQoldy = dfxdiQy(it);
    dfydrQoldx = dfydrQx(it);
    dfydiQoldx = dfydiQx(it);
    
end

contactflag = 1;

Fx = ffttho(ftx(length(ftx)/2+1:end),1,'nfft',length(ftx)/2); %only use second period for FFT
Fy = ffttho(fty(length(fty)/2+1:end),1,'nfft',length(fty)/2); %only use second period for FFT
dFxdrQx = ffttho(dfxdrQx(length(dfxdrQx)/2+1:end),1,'nfft',length(dfxdrQx)/2); %only use second period for FFT
dFxdiQx = ffttho(dfxdiQx(length(dfxdiQx)/2+1:end),1,'nfft',length(dfxdiQx)/2); %only use second period for FFT
dFydrQy = ffttho(dfydrQy(length(dfydrQy)/2+1:end),1,'nfft',length(dfydrQy)/2); %only use second period for FFT
dFydiQy = ffttho(dfydiQy(length(dfydiQy)/2+1:end),1,'nfft',length(dfydiQy)/2); %only use second period for FFT
dFydrQx = ffttho(dfydrQx(length(dfydrQx)/2+1:end),1,'nfft',length(dfydrQx)/2); %only use second period for FFT
dFydiQx = ffttho(dfydiQx(length(dfydiQx)/2+1:end),1,'nfft',length(dfydiQx)/2); %only use second period for FFT
dFxdrQy = ffttho(dfxdrQy(length(dfxdrQy)/2+1:end),1,'nfft',length(dfxdrQy)/2); %only use second period for FFT
dFxdiQy = ffttho(dfxdiQy(length(dfxdiQy)/2+1:end),1,'nfft',length(dfxdiQy)/2); %only use second period for FFT

Foutx  = [2*real(Fx(length(Fx)/2+2)); -2*imag(Fx(length(Fx)/2+2))]; %first harmonic
Fouty  = [2*real(Fy(length(Fy)/2+2)); -2*imag(Fy(length(Fy)/2+2))]; %first harmonic
dFxdQoutx = [2*real(dFxdrQx(length(dFxdrQx)/2+2)), 2*real(dFxdiQx(length(dFxdiQx)/2+2)); -2*imag(dFxdrQx(length(dFxdrQx)/2+2)), -2*imag(dFxdiQx(length(dFxdiQx)/2+2))];
dFydQouty = [2*real(dFydrQy(length(dFydrQy)/2+2)), 2*real(dFydiQy(length(dFydiQy)/2+2)); -2*imag(dFydrQy(length(dFydrQy)/2+2)) , -2*imag(dFydiQy(length(dFydiQy)/2+2))];
dFxdQouty = [2*real(dFxdrQy(length(dFxdrQy)/2+2)), 2*real(dFxdiQy(length(dFxdiQy)/2+2)); -2*imag(dFxdrQy(length(dFxdrQy)/2+2)), -2*imag(dFxdiQy(length(dFxdiQy)/2+2))];
dFydQoutx = [2*real(dFydrQx(length(dFydrQx)/2+2)), 2*real(dFydiQx(length(dFydiQx)/2+2)); -2*imag(dFydrQx(length(dFydrQx)/2+2)), -2*imag(dFydiQx(length(dFydiQx)/2+2))];

fnl = [Foutx;Fouty;0;0;0;0];
R = S*qamp - fl - fnl;
dFdQ = kron(diag([1,0,0,0]),dFxdQoutx)+kron(diag([0,1,0,0]),dFydQouty);
dFdQ(1:2,3:4)=dFxdQouty;
dFdQ(3:4,1:2)=dFydQoutx;
J = S-dFdQ;


