% 2 DOF simple friction damped blade model. DOF 1 is disk with friction
% damping, DOF 2 is blade with excitation.
%==========================================================================

clear
clc
%load('/home/SERVER-hoffmann/sync/03_Matlab/AFT/nonlinearfrictionforce/hysteresis.mat')

%% UI
%frequency
freqmin = 0;
freqmax = 200;
freqstep = 1;

% System definition q = [q_d,c; q_d,s; q_b,c; q_b,s]
md = 1;         %disk massen
mb = 0.3;       %blade massen
kgdx = 100e3;  %ground- disk -x
kgdy = 100e3;  %ground- disk -y
kgbx = 0;  %ground- blade -x
kgby = 0;  %ground- blade -y
kdb = 100e3;  %disk- blade
ddx = 5;   %structural damping  disk-x
ddy = 5;   %structural damping  disk-y
dbx = 5;   %structural damping  blade-x
dby = 5;   %structural damping  blade-y

Msingle = diag([md,mb]);        % ungekoppelte massenmatrix
M = kron(Msingle,eye(4));       % 4x4
%Ksingle = [kgd+kdb, -kdb; -kdb, kdb];
%K = kron(Ksingle,eye(2));
Ksingle = diag([kgdx+kdb,kgdx+kdb,kgdy+kdb,kgdy+kdb,kdb,kdb,kgby+kdb,kgby+kdb]);
Keck = kron([0,1;1,0],diag([-kdb,-kdb,-kdb,-kdb]));
K = Ksingle + Keck;


Dsingle = zeros(8,8);
Dsingle(1,2) = ddx;
Dsingle(2,1) = -ddx;
Dsingle(3,4) = ddy;
Dsingle(4,3) = -ddy;
Dsingle(5,6) = dbx;
Dsingle(6,5) = -dbx;
Dsingle(7,8) = dby;
Dsingle(8,7) = -dby;
D = Dsingle;
%D = kron(eye(2),Dsingle);
% Dsingle = [d,0;0,d];
% D = kron(Dsingle,eye(2));
% [V E] = eig(Ksingle,Msingle);
% E = sqrt(E)/2/pi;

% Numerical parameters




% Contact Parameters
kt = 50e6;
fnormal = 8;
mu = 1;
mu_x = 1;
mu_y = 1;

%% Forcing

fl  = [0; 0; 0; 0; 0; 10; 0; 10]; % f = [......; f_b_x,c; f_b_x,s; f_b_y,c; f_b_y,s]
fnl = [0; 0; 0; 0; 0; 0; 0; 0;];

%% Main Loop
freqrange = freqmin:freqstep:freqmax;
Q0 = [0; 0; 0; 0; 0; 0; 0; 0];  % q_disk_x_cos; q_disk_x_sin; q_disk_y_cos; q_disk_y_sin
                                % q_blade_x_cos; q_blade_x_sin; q_blade_y_cos; q_blade_y_sin
Qvec = zeros(length(freqrange),length(Q0)/2);
status = zeros(length(freqrange),1);
tic
% -------------------------------------------------------------------------    
% frequency loop    
% -------------------------------------------------------------------------
    Soptfsolve = optimoptions(@fsolve,...
    'Display','none', ...
    'TolX',1.0E-3, ...
    'Tolfun',1.0E-3, ...
    'MaxIter',100,...
    'Diagnostics','off', ...
    'FunValCheck','off',...
    'DerivativeCheck','off',...
    'Jacobian','on');
for ifreq = 1: length(freqrange) % frequencyloop
    freq = freqrange(ifreq);
    S = calc_dynamic_stiffness(freq,K,M,D); %freq ==0 S==K

    debug = 0;
    [Q,R,ikeyexit,Soutput,J] = fsolve(@(Q) func_residualaft2mal1D(Q,S,fl,kt,fnormal,mu_x,mu_y,debug),Q0,Soptfsolve);
%     [Q,R,ikeyexit,Soutput,J] = fsolve(@(Q) func_residualaft2D(Q,S,fl,kt,fnormal,mu,debug),Q0,Soptfsolve);
    fprintf([' ikeyexit = ',num2str(ikeyexit),'\n']);
    res = norm(R);
    init = Soutput.iterations;
    infeval = Soutput.funcCount;
%     disp(Soutput.message);
    % Show user information
    fprintf('Number of iterations               : %12i [-]\n',init);
    if ikeyexit < 1
        fprintf('no convergence');
    end
    tmaxQdx = 1/2/pi/freq*atan2(Q(2),Q(1));
    tmaxQdy = 1/2/pi/freq*atan2(Q(4),Q(3));% find amplitude from 
    tmaxQbx = 1/2/pi/freq*atan2(Q(6),Q(5));
    tmaxQby = 1/2/pi/freq*atan2(Q(8),Q(7));
    Qvec_disk(ifreq,:) = [Q(1)*cos(2*pi*freq*tmaxQdx)+Q(2)*sin(2*pi*freq*tmaxQdx), Q(3)*cos(2*pi*freq*tmaxQdy)+Q(4)*sin(2*pi*freq*tmaxQdy)];
    Qvec_blade(ifreq,:) = [Q(5)*cos(2*pi*freq*tmaxQbx)+Q(6)*sin(2*pi*freq*tmaxQbx), Q(7)*cos(2*pi*freq*tmaxQby)+Q(8)*sin(2*pi*freq*tmaxQby)];
    
    Q0 = Q;
end
% -------------------------------------------------------------------- *freq*tmaxF);


% -------------------------------------------------------------------------    
% end frequency loop    
% -------------------------------------------------------------------------
t = toc;
%% Postprocessing


figure(105)
hold on
semilogy(freqrange,abs(Qvec_disk));
str101 = ['Maximal Displacement of Disk with \mu=', num2str(mu), ' and Fn=', num2str(fnormal)];
title(str101);
legend('Disk\_x','Disk\_y');
xlabel('Frequency(Hz)')
ylabel('Displacement(m)')
grid on
figure(106)
hold on
semilogy(freqrange,abs(Qvec_blade));
str102 = ['Maximal Displacement of Blade with \mu=', num2str(mu), ' and Fn=', num2str(fnormal)];
title(str102);
legend('Blade\_x','Blade\_y');
xlabel('Frequency(Hz)')
ylabel('Displacement(m)')
grid on
%ylabel('Displacement [10^{-6} m]')
%xlabel('Frequency [Hz]')
%legend(legtext,'Fontsize',8);
%set(gcf,'Color','w','Units','centimeter','Position',[5,5,15,10])
%set(gca,'Fontsize',8)


%%