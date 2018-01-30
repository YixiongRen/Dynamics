clear;
clc;
clf;
%% initializing

miu = 0.3;  %      empfindlich miu Inversely proportional to u
q0 = 100;           % Max_Ff
q2 = q0/2;           % Max_Ff    q Inversely proportional to u
E = 21000000000;        % beeinfluss die Verschiebung  E Inversely proportional to u
A = 0.0001;         % beeinfluss die Verschiebung  A Inversely proportional to u
L = 0.002;   %      empfindlich L proportional to u
F_amp = 50;         
Omega = 30;        
%% Loading force
% delta = [];
% for t = 0:0.001:pi/2
%     F_t = F_amp*sin(t);
%     P_F = [-4*miu*q2/3/L^2,2*miu*q2/L,miu*q0,-F_t];
%     delta(end+1) = fzero(@(x) -4*miu*q2/3/L^2*x^3+2*miu*q2/L*x^2+miu*q0*x-F_t,-1);
% %     Roots = roots(P_F);
% %     x(end+1) = Roots(1);
% end
% delta(:) = delta(:)-delta(1);
% u = [];
% for t = 1:length(delta)
%     u(end+1) = miu*delta(t)^2/E/A*(q0/2+q2*(4*delta(t)/3/L-delta(t)^2/L^2));
% end
% t = 0:0.001:pi/2;
% ft = F_amp*sin(t);
% plot(u,ft,'color','r');
% u_amp = u(end);
% F_amp = ft(end);
% delta_amp = delta(end);
% fun = inline('-4*miu*q2/3/L^2*x^3+2*miu*q2/L*x^2+miu*q0*x-F_amp','miu','q2','L','q0','F_amp','x');
% fun=inline('45*x+30/0.25*(1.5*x^2-2*x^3)-50','x');  % also replace the F= 50 here
delta_amp = fzero(@Func,[-10,10]);



u_amp = miu*delta_amp^2/E/A*(q0/2+q2*(4*delta_amp/3/L-delta_amp^2/L^2));
delta = 0:-0.0001:delta_amp;
u = [];
F = [];
for i = 1:length(delta)
    u(end+1) = miu*delta(i)^2/E/A*(q0/2+q2*(4/3/L*delta(i)-1/L^2*delta(i)^2));
   
end

for i = 1:length(delta)
    F(end+1) = miu*q0*delta(i)+2*miu*q2/L*delta(i)^2-4*miu*q2/3/L^2*delta(i)^3;
   
end
plot(u,F,'color','r','LineWidth',2);
%% F decreasing
delta = 0:-0.0001:delta_amp;
u_d = [];
F_d = [];
for i = 1:length(delta)
    u_d(end+1) = miu/E/A*(q0*(delta_amp^2-2*delta(i)^2)/2+q2*(4/3/L*(delta_amp^3-2*delta(i)^3)-(delta_amp^4-2*delta(i)^4)/L^2));
      
end

for i = 1:length(delta)
    F_d(end+1) = miu*q0*(delta_amp-2*delta(i))+2*miu*q2/3/L^2*(3*L*delta_amp^2-2*delta_amp^3+4*delta(i)^3-6*L*delta(i)^2);
   
end
hold on;
plot(u_d,F_d,'color','b','LineWidth',2);

%% F increasing

delta = 0:-0.0001:delta_amp;
u_i = [];
F_i = [];
for i = 1:length(delta)
    u_i(end+1) = miu/E/A*(q0*(2*delta(i)^2-delta_amp^2)/2+q2*(4/3/L*(2*delta(i)^3-delta_amp^3)-(2*delta(i)^4-delta_amp^4)/L^2));
   
end

for i = 1:length(delta)
    F_i(end+1) = miu*q0*(2*delta(i)-delta_amp)+2*miu*q2/3/L^2*(-3*L*delta_amp^2+2*delta_amp^3-4*delta(i)^3+6*L*delta(i)^2);
   
end
hold on;
plot(u_i,F_i,'color','c','LineWidth',2);

xlabel('u');
ylabel('F');



