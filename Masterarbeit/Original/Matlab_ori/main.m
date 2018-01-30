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
kgd = 100e3;  %ground- disk
kdb = 100e3;  %disk- blade
d = 5;   %structural damping

Msingle = diag([md,mb]);        % ungekoppelte massenmatrix
M = kron(Msingle,eye(2));       % 4x4
Ksingle = [kgd+kdb, -kdb; -kdb, kdb];
K = kron(Ksingle,eye(2));

Dsingle = [0, d; -d, 0];
D = kron(eye(2),Dsingle);
% Dsingle = [d,0;0,d];
% D = kron(Dsingle,eye(2));
% [V E] = eig(Ksingle,Msingle);
% E = sqrt(E)/2/pi;

% Numerical parameters

epsilon = 1e-8;
iitmax = 100;
newdampmax = 4;

% Contact Parameters
kt = 50e6;
fnormal = 35;
mu = 0.5;

%% Forcing

fl  = [0; 0; 20; 0]; % f = [f_d,c; f_d,s; f_b,c; f_b,s]
fnl = [0; 0; 0; 0];

%% Main Loop
freqrange = freqmin:freqstep:freqmax;
Q = [1; 0; 1; 0];  % q_disk_cos; q_disk_sin; q_blade_cos; q_blade_sin
Qvec = zeros(length(freqrange),length(Q)/2);
status = zeros(length(freqrange),1);
tic
% -------------------------------------------------------------------------    
% frequency loop    
% -------------------------------------------------------------------------
for ifreq = 1: length(freqrange) % frequencyloop
    freq = freqrange(ifreq);
    S = calc_dynamic_stiffness(freq,K,M,D);
% -------------------------------------------------------------------------    
% iteration loop    
% -------------------------------------------------------------------------
    for iit = 1: iitmax 
        convergence = 0;
        if iit == iitmax
            error('Maximum number of iterations reached');
        end
        
        if iit == 1
            newdamp = 0;
        else
            newdamp = newdampmax;
        end
        
        %calculate harmonics of nonlinear friction force
        Qcont = [Q(1),Q(2)];
        [fnlcont,dFdQ, contactflag, contactstatus] = calc_nonlinfricforce (Qcont,kt,fnormal,mu,0);
        fnl = [fnlcont;0;0];
        %updating dynamic stiffness matrix
        
        S = calc_dynamic_stiffness(freq,K,M,D);
        
        % calculate residual
        R = S*Q - fl - fnl;
        Qold = Q;
        normrso = 0;
% -------------------------------------------------------------------------    
% newton iteration loop    
% -------------------------------------------------------------------------
        % calculate newton step
        J = S-kron(diag([1,0]),dFdQ);
        
        dQ = -J\R; % minus as in wikipedia. S is here an approximation for jacobimatrix.
        
        for inew = 0:newdamp
            
            Q = Qold + dQ/2^inew;
            
            norm.Q = abs(Q);
            norm.dQ = abs(dQ)/2^inew;
            
            % Check the convergence criterion and proceed to next frequency step
            if (norm.dQ/norm.Q) < epsilon
                format1 = ['\n Relative change in amplitude satisfies tolerance ' ...
                    'after %i equilibrium iterations !\n'];
                fprintf(1, format1, iit);
                convergence = 1;
                break;
            end
            %calculate harmonics of nonlinear friction force
            
            Qcont = [Q(1),Q(2)];
            [fnlcont,dFdQ, contactflag, contactstatus] = calc_nonlinfricforce (Qcont,kt,fnormal,mu,0);
            fnl = [fnlcont;0;0];
        
            %updating dynamic stiffness matrix

            S = calc_dynamic_stiffness(freq,K,M,D);
        
            % calculate residual
            R = S*Q - fl - fnl;
            norm.res = R'*R;
            
            if ((norm.res > normrso) && (inew > 0))
                %Accept the damped Newton step and advance to the next iteration step
                Q = Qold + dQ/(2^(inew-1));
                break;
            else
                normrso = norm.res;
            end
            
            
        end
% -------------------------------------------------------------------------    
% end newton iteration loop    
% -------------------------------------------------------------------------
        if convergence == 1
            format2 = ['\n Convergence accomplished for frequency ',...
                       '%.4f Hz!\n'];
            fprintf(1, format2, freq);
            break
        end
    end 
% -------------------------------------------------------------------------    
% end iteration loop    
% -------------------------------------------------------------------------
    tmaxQd = 1/2/pi/freq*atan2(Q(2),Q(1));
    tmaxQb = 1/2/pi/freq*atan2(Q(4),Q(3));% find amplitude from 
    Qvec(ifreq,:) = [Q(1)*cos(2*pi*freq*tmaxQd)+Q(2)*sin(2*pi*freq*tmaxQd), Q(3)*cos(2*pi*freq*tmaxQb)+Q(4)*sin(2*pi*freq*tmaxQb)];
%     tmaxF = 1/2/pi/freq*atan2(fnl(2),fnl(1));% find amplitude from 
%     Fvec(ifreq) = fnl(1)*cos(2*pi*freq*tmaxF)+fnl(2)*sin(2*pi*freq*tmaxF);

end 
% -------------------------------------------------------------------------    
% end frequency loop    
% -------------------------------------------------------------------------
t = toc;
%% Postprocessing


figure(101)
semilogy(freqrange,abs(Qvec)); hold on
%ylabel('Displacement [10^{-6} m]')
%xlabel('Frequency [Hz]')
%legend(legtext,'Fontsize',8);
%set(gcf,'Color','w','Units','centimeter','Position',[5,5,15,10])
%set(gca,'Fontsize',8)
hold on

%%