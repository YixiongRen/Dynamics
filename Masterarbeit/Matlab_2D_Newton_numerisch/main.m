% 2 DOF simple friction damped blade model. DOF 1 is disk with friction
% damping, DOF 2 is blade with excitation.
%==========================================================================


%load('/home/SERVER-hoffmann/sync/03_Matlab/AFT/nonlinearfrictionforce/hysteresis.mat')

%% UI
%frequency
freqmin = 0;
freqmax = 200;
freqstep = 5;

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

Msingle = diag([md,md,mb,mb]);        % ungekoppelte massenmatrix
M = kron(Msingle,eye(2));       % 4x4
%Ksingle = [kgd+kdb, -kdb; -kdb, kdb];
%K = kron(Ksingle,eye(2));
% Ksingle = diag([kgdx+kdb,kgdx+kdb,kgdy+kdb,kgdy+kdb,kdb,kdb,kgby+kdb,kgby+kdb]);
% Keck = kron([0,1;1,0],diag([-kdb,-kdb,-kdb,-kdb]));
Ksingle = diag([kgdx+kdb,kgdy+kdb,kdb,kgby+kdb]);
Ksingle(1,3) = -kdb;
Ksingle(2,4) = -kdb;
Ksingle(3,1) = -kdb;
Ksingle(4,2) = -kdb;
K = kron(Ksingle,eye(2));
% K = Ksingle + Keck;


% Dsingle = zeros(8,8);
% Dsingle(1,2) = ddx;
% Dsingle(2,1) = -ddx;
% Dsingle(3,4) = ddy;
% Dsingle(4,3) = -ddy;
% Dsingle(5,6) = dbx;
% Dsingle(6,5) = -dbx;
% Dsingle(7,8) = dby;
% Dsingle(8,7) = -dby;
% D = Dsingle;
Dsingle = [0, ddx,0,0; -ddy, 0,0,0; 0,0,0,dbx; 0,0,-dby,0];
D = kron(eye(2),Dsingle);
%D = kron(eye(2),Dsingle);
% Dsingle = [d,0;0,d];
% D = kron(Dsingle,eye(2));
% [V E] = eig(Ksingle,Msingle);
% E = sqrt(E)/2/pi;

% Numerical parameters

epsilon = 1e-8;
iitmax = 1000;
newdampmax = 40;

% Contact Parameters
kt = 50e6;
fnormal = 15;
mu = 1;

%% Forcing

fl  = [0; 0; 0; 0; 0; 0; 0; 10]; % f = [......; f_b_x,c; f_b_x,s; f_b_y,c; f_b_y,s]
fnl = [0; 0; 0; 0; 0; 0; 0; 0;];

%% Main Loop
freqrange = freqmin:freqstep:freqmax;
contactstatu = zeros(length(freqrange),2);
Q = [0; 0; 0; 0; 0; 0; 0; 0];  % q_disk_x_cos; q_disk_x_sin; q_disk_y_cos; q_disk_y_sin
                                % q_blade_x_cos; q_blade_x_sin; q_blade_y_cos; q_blade_y_sin
Qvec = zeros(length(freqrange),length(Q)/2);
status = zeros(length(freqrange),1);
tic
% -------------------------------------------------------------------------    
% frequency loop    
% -------------------------------------------------------------------------
for ifreq = 1: length(freqrange) % frequencyloop
    freq = freqrange(ifreq);
    S = calc_dynamic_stiffness(freq,K,M,D); %freq ==0 S==K
% -------------------------------------------------------------------------    
% iteration loop    
% -------------------------------------------------------------------------
    for iit = 1: iitmax 
        convergence = 0;
        if iit == iitmax
            disp('Maximum number of iterations reached');
        end
        
        if iit == 1
            newdamp = 0;
        else
            newdamp = newdampmax;
        end
        
        %calculate harmonics of nonlinear friction force
        Qcont = [Q(1),Q(2),Q(3),Q(4)];  % disk_x_cos, disk_x_sin, disk_y_cos, disk_y_sin
        %[fnlcontx,fnlconty,dFdQx,dFdQy, contactflag, contactstatus] = calc_nonlinfricforce (Qcont,kt,fnormal,mu,0,0);
        %contactstatu(ifreq,end+1)=contactstatus;
        [Fout] = calc_nlforce(Qcont,kt,fnormal,mu);
        [jac,err] = jacobianest(@(qamp) calc_nlforce(qamp,kt,fnormal,mu),Qcont);
        
        
        fnl = [Fout(1:2)';Fout(3:4)';0;0;0;0];
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
        dFxdQx = jac(1:2,1:2);
        dFydQy = jac(3:4,3:4);
        dFxdQy = jac(1:2,3:4);
        dFydQx = jac(3:4,1:2);
        dFdQ = kron(diag([1,0,0,0]),dFxdQx)+kron(diag([0,1,0,0]),dFydQy);
        dFdQ(1:2,3:4) = dFxdQy;
        dFdQ(3:4,1:2) = dFydQx;
        J = S-dFdQ;
        
%        dQ = -pinv(J)*R; 
        dQ = -J\R;               % minus as in wikipedia. S is here an approximation for jacobimatrix.
        
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
                %[fnlcontx,fnlconty,dFdQx,dFdQy, contactflag, contactstatus] = calc_nonlinfricforce (Qcont,kt,fnormal,mu,0,1);
                [Fout] = calc_nlforce(Qcont,kt,fnormal,mu);
                [jac,err] = jacobianest(@(qamp) calc_nlforce(qamp,kt,fnormal,mu),Qcont);
                dFxdQx = jac(1:2,1:2);
                dFydQy = jac(3:4,3:4);
                dFxdQy = jac(1:2,3:4);
                dFydQx = jac(3:4,1:2);
                break;
            end
            %calculate harmonics of nonlinear friction force
            
            Qcont = [Q(1),Q(2),Q(3),Q(4)];
            %[fnlcontx,fnlconty,dFdQx,dFdQy, contactflag, contactstatus] = calc_nonlinfricforce (Qcont,kt,fnormal,mu,0,0);
            %contactstatu(ifreq,end+1)=contactstatus;
            [Fout] = calc_nlforce(Qcont,kt,fnormal,mu);
            [jac,err] = jacobianest(@(qamp) calc_nlforce(qamp,kt,fnormal,mu),Qcont);
            dFxdQx = jac(1:2,1:2);
            dFydQy = jac(3:4,3:4);
            dFxdQy = jac(1:2,3:4);
            dFydQx = jac(3:4,1:2);
            fnl = [Fout(1:2)';Fout(3:4)';0;0;0;0];
        
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
    tmaxQdx = 1/2/pi/freq*atan2(Q(2),Q(1));
    tmaxQdy = 1/2/pi/freq*atan2(Q(4),Q(3));% find amplitude from 
    tmaxQbx = 1/2/pi/freq*atan2(Q(6),Q(5));
    tmaxQby = 1/2/pi/freq*atan2(Q(8),Q(7));
    Qvec_disk(ifreq,:) = [Q(1)*cos(2*pi*freq*tmaxQdx)+Q(2)*sin(2*pi*freq*tmaxQdx), Q(3)*cos(2*pi*freq*tmaxQdy)+Q(4)*sin(2*pi*freq*tmaxQdy)];
    Qvec_blade(ifreq,:) = [Q(5)*cos(2*pi*freq*tmaxQbx)+Q(6)*sin(2*pi*freq*tmaxQbx), Q(7)*cos(2*pi*freq*tmaxQby)+Q(8)*sin(2*pi*freq*tmaxQby)];
%     tmaxF = 1/2/pi/freq*atan2(fnl(2),fnl(1));% find amplitude from 
%     Fvec(ifreq) = fnl(1)*cos(2*pi*freq*tmaxF)+fnl(2)*sin(2*pi*freq*tmaxF);

end 
% -------------------------------------------------------------------------    
% end frequency loop    
% -------------------------------------------------------------------------
t = toc;
%% Postprocessing


figure(101)
hold on
semilogy(freqrange,abs(Qvec_disk));
str101 = ['Maximal Displacement of Disk with \mu=', num2str(mu), ' and Fn=', num2str(fnormal)];
title(str101);
legend('Disk\_x','Disk\_y');
xlabel('Frequency(Hz)')
ylabel('Displacement(m)')
grid on
figure(102)
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
% figure(103)
% imagesc(contactstatu)
% colorbar
% disp(length(find(contactstatu==2)));
%%