function [Fout, dFdQout, contactflag, contactstatus] = calc_nonlinfricforce (qamp,kt,fn,mu,debug)
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
q = qamp(1)*cos(qt)+qamp(2)*sin(qt);
N = length(q);

% Solving process

%initialization
y = zeros(length(qt),1);  % old Q
ft = zeros(length(qt),1);
dfdrQold = 0;   % real part
dfdiQold = 0;   % imaginary part
dfdrQ = zeros(length(qt),1);
dfdiQ = zeros(length(qt),1);

contactstatus = 0;
%========================================================================== 
% main integration loop
%========================================================================== 
for it = 2: length(qt)
% ---------------Slip in positive direction--------------------------------
    if kt*(q(it) - y(it-1))- mu*fn > 0 % slip pos
        ft(it) = -mu*fn;
        dfdrQ(it) = 0;
        dfdiQ(it) = 0;
        y(it) = q(it) + ft(it)/kt;
        contactstatus = 1;
% ---------------Slip in negative direction--------------------------------
    elseif -kt*(q(it) - y(it-1))- mu*fn > 0 % slip neg
        ft(it) = +mu*fn;
        dfdrQ(it) = 0;
        dfdiQ(it) = 0;
        y(it) = q(it) + ft(it)/kt;
        contactstatus = 1;
% ---------------General stick---------------------------------------------
    elseif kt*(q(it) - y(it-1))- mu*fn <= 0
        ft(it) = -kt*(q(it) - y(it-1));
        dfdrQ(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQold; % the factor 2*2*pi is hard coded due to 2 periods!
        dfdiQ(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQold;
        y(it) = y(it-1);
    end
    dfdrQold = dfdrQ(it);
    dfdiQold = dfdiQ(it);

    
end

contactflag = 1;

F = ffttho(ft(length(ft)/2+1:end),1,'nfft',length(ft)/2); %only use second period for FFT
dFdrQ = ffttho(dfdrQ(length(dfdrQ)/2+1:end),1,'nfft',length(dfdrQ)/2); %only use second period for FFT
dFdiQ = ffttho(dfdiQ(length(dfdiQ)/2+1:end),1,'nfft',length(dfdiQ)/2); %only use second period for FF


Fout  = [2*real(F(length(F)/2+2)); -2*imag(F(length(F)/2+2))]; %first harmonic
dFdQout = [2*real(dFdrQ(length(dFdrQ)/2+2)), 2*real(dFdiQ(length(dFdiQ)/2+2)); -2*imag(dFdrQ(length(dFdrQ)/2+2)), -2*imag(dFdiQ(length(dFdiQ)/2+2))];



%% Debugging section
if debug == 1
%     figure(200)
%     stem(abs(2*F(length(F)/2+1:end)))
    
%     figure(201)
%     plot(qt,q,qt,y)
%     legend('Excitation','Response')
    N = 2048;
    Funit = [zeros(N/2-1,1);F(1024:1:1026);zeros(N/2-2,1)];  
    
    f = real(ffttho(Funit,-1));

    figure(202)
    plot(qt(1:2048),q(2049:end),qt(1:2048),ft(2049:end)',qt(1:2048),f);
    legend({'Verschiebung','Kraft','Kraft monoharmonisch'})
    
%      figure(203)
%      plot(qt,dfdrQ,qt,dfdiQ);
%      legend('Re(dF/Re(dQ))','Im(dF/Re(dQ))')
%    
    figure(204)
    plot(q(2049:end)',ft(2049:end),q(2049:end),f);
    legend({'Voll','Monoharmonisch'});
end

