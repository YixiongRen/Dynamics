function [R,J] = func_residualaft(qamp, S, fl, kt, fn, mu, debug)
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
N = length(qx);

% Solving process

%initialization
x = zeros(length(qt),1);  % old Q
y = zeros(length(qt),1);  % old Q
ftx = zeros(length(qt),1);
fty = zeros(length(qt),1);
dfdrQoldx = 0;   % real part
dfdiQoldx = 0;   % imaginary part
dfdrQoldy = 0;   % real part
dfdiQoldy = 0;   % imaginary part
dfdrQx = zeros(length(qt),1);
dfdiQx = zeros(length(qt),1);
dfdrQy = zeros(length(qt),1);
dfdiQy = zeros(length(qt),1);

contactstatus = 0;
%========================================================================== 
% main iteration loop
%========================================================================== 
for it = 2: length(qt)
% ---------------Slip in all positive directions--------------------------------
    if kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn > 0 && qx(it) > x(it-1) && qy(it) > y(it-1) % slip pos
        ftx(it) = -mu*fn;
        fty(it) = -mu*fn;
        dfdrQx(it) = 0;
        dfdiQx(it) = 0;
        dfdrQy(it) = 0;
        dfdiQy(it) = 0;
        x(it) = qx(it) - abs(ftx(it)/kt*cos(atan((qy(it)-y(it-1))/(qx(it))-x(it-1)))); % x == qx- feder-extension
        y(it) = qy(it) - abs(fty(it)/kt*sin(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        contactstatus = 1;
% ---------------Slip in negative y direction--------------------------------
    elseif kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn > 0 && qx(it) > x(it-1) && qy(it) <= y(it-1) % slip neg
        ftx(it) = -mu*fn;
        fty(it) = +mu*fn;
        dfdrQx(it) = 0;
        dfdiQx(it) = 0;
        dfdrQy(it) = 0;
        dfdiQy(it) = 0;
        x(it) = qx(it) - abs(ftx(it)/kt*cos(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        y(it) = qy(it) + abs(fty(it)/kt*sin(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        contactstatus = 2;
% ---------------Slip in negative x direction--------------------------------
    elseif kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn > 0 && qx(it) <= x(it-1) && qy(it) > y(it-1) % slip neg
        ftx(it) = +mu*fn;
        fty(it) = -mu*fn;
        dfdrQx(it) = 0;
        dfdiQx(it) = 0;
        dfdrQy(it) = 0;
        dfdiQy(it) = 0;
        x(it) = qx(it) + abs(ftx(it)/kt*cos(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        y(it) = qy(it) - abs(fty(it)/kt*sin(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        contactstatus = 3;
% ---------------Slip in all negative directions--------------------------------
    elseif kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn > 0 && qx(it) <= x(it-1) && qy(it) <= y(it-1) % slip pos
        ftx(it) = +mu*fn;
        fty(it) = +mu*fn;
        dfdrQx(it) = 0;
        dfdiQx(it) = 0;
        dfdrQy(it) = 0;
        dfdiQy(it) = 0;
        x(it) = qx(it) + abs(ftx(it)/kt*cos(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        y(it) = qy(it) + abs(fty(it)/kt*sin(atan((qy(it)-y(it-1))/(qx(it))-x(it-1))));
        contactstatus = 4;
% ---------------General stick---------------------------------------------
    elseif kt*sqrt((qx(it) - x(it-1))^2 + (qy(it)-y(it-1))^2) - mu*fn <= 0
        ftx(it) = -kt*(qx(it) - x(it-1));
        fty(it) = -kt*(qy(it) - y(it-1));
        dfdrQx(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldx; % the factor 2*2*pi is hard coded due to 2 periods!
        dfdiQx(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldx;
        dfdrQy(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldy; % the factor 2*2*pi is hard coded due to 2 periods!
        dfdiQy(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldy;
        x(it) = x(it-1);
        y(it) = y(it-1);
    end
    dfdrQoldx = dfdrQx(it);
    dfdiQoldx = dfdiQx(it);
    dfdrQoldy = dfdrQy(it);
    dfdiQoldy = dfdiQy(it);

    
end

contactflag = 1;

Fx = ffttho(ftx(length(ftx)/2+1:end),1,'nfft',length(ftx)/2); %only use second period for FFT
Fy = ffttho(fty(length(fty)/2+1:end),1,'nfft',length(fty)/2); %only use second period for FFT
dFdrQx = ffttho(dfdrQx(length(dfdrQx)/2+1:end),1,'nfft',length(dfdrQx)/2); %only use second period for FFT
dFdiQx = ffttho(dfdiQx(length(dfdiQx)/2+1:end),1,'nfft',length(dfdiQx)/2); %only use second period for FF
dFdrQy = ffttho(dfdrQy(length(dfdrQy)/2+1:end),1,'nfft',length(dfdrQy)/2); %only use second period for FFT
dFdiQy = ffttho(dfdiQy(length(dfdiQy)/2+1:end),1,'nfft',length(dfdiQy)/2); %only use second period for FF

Foutx  = [2*real(Fx(length(Fx)/2+2)); 2*imag(Fx(length(Fx)/2+2))]; %first harmonic
Fouty  = [2*real(Fy(length(Fy)/2+2)); 2*imag(Fy(length(Fy)/2+2))]; %first harmonic
dFdQoutx = [2*real(dFdrQx(length(dFdrQx)/2+2)), 2*real(dFdiQx(length(dFdiQx)/2+2)); 2*imag(dFdrQx(length(dFdrQx)/2+2)), 2*imag(dFdiQx(length(dFdiQx)/2+2))];
dFdQouty = [2*real(dFdrQy(length(dFdrQy)/2+2)), 2*real(dFdiQy(length(dFdiQy)/2+2)); 2*imag(dFdrQy(length(dFdrQy)/2+2)), 2*imag(dFdiQy(length(dFdiQy)/2+2))];

fnl = [Foutx;Fouty;0;0;0;0];
R = S*qamp - fl - fnl;
dFdQ = kron(diag([1,0,0,0]),dFdQoutx)+kron(diag([0,1,0,0]),dFdQouty);
J = S-dFdQ;






%% Debugging section
if debug == 1
%     figure(200)
%     stem(abs(2*F(length(F)/2+1:end)))
    
%     figure(201)
%     plot(qt,q,qt,y)
%     legend('Excitation','Response')
    N = 2048;
    Funitx = [zeros(N/2-1,1);Fx(1024:1:1026);zeros(N/2-2,1)];  
    
    fx = real(ffttho(Funitx,-1));

    figure(202)
    plot(qt(1:2048),qx(2049:end),qt(1:2048),ftx(2049:end)',qt(1:2048),fx);
    legend({'Verschiebung','Kraft','Kraft monoharmonisch'})
    
%      figure(203)
%      plot(qt,dfdrQ,qt,dfdiQ);
%      legend('Re(dF/Re(dQ))','Im(dF/Re(dQ))')
%    
    figure(204)
    plot(qx(2049:end)',ftx(2049:end),qx(2049:end),fx);
    legend({'Voll','Monoharmonisch'});
end