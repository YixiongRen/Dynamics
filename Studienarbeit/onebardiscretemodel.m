clear;
clc;
clf;
%% initializing

miu = 0.15;  %
q0 = 300;
q2 = 300;
E = 2100000;
A = 0.0003;
L = 0.5;   %
F_amp = 50;
Omega = 300;
N = 500;
delta_x = L/N;
syms x;
Sig = [];
u_summe = [];
Ff_summe = [];
for k = 1:N
    Sig(end+1) = 0;
end

%% 
for t = 0:0.01:pi*5/2
    F_sum = 0;              % generate empty arrays
    delta_Ff = [];
    delta_Ypsilon = [];
    Ypsilon = [];
    u = [];
    F = F_amp*sin(t);

    
    if F == 0
      
    else
        for i = 0:(N-1)
            if F > 0 && i==0        % first node
                       
                F_sum = F_sum + double(int( q0+q2*4/L^2*(x*L-x^2), delta_x*i, delta_x*(i+1) ));
                Sig(i+1) = 1;
            elseif F < 0 && i==0    % first node
        
                F_sum = F_sum - double(int( q0+q2*4/L^2*(x*L-x^2), delta_x*i, delta_x*(i+1) ));
                Sig(i+1) = -1;
            else                    % the rest
                if F_sum >= F && F > 0      % Ff is enough large
                    break;
                elseif F_sum <= F && F < 0  % Ff is enough small
                    break;
                elseif F_sum <= F && F > 0  % Ff is not enough large
                    F_sum = F_sum + double(int( q0+q2*4/L^2*(x*L-x^2), delta_x*i, delta_x*(i+1) ));
                    Sig(i+1) = 1;
                elseif F_sum >= F && F < 0  % Ff is not enough small
                    F_sum = F_sum - double(int( q0+q2*4/L^2*(x*L-x^2), delta_x*i, delta_x*(i+1) ));
                    Sig(i+1) = -1;  
                else    
                end
            end         
        
            
        end  % finish generated Sig(:) for a unique F.
        
        
        % should here calculate the delta_Ff and delta_Ypsilon from Sig(:)
        for j = 1:length(Sig)
            if Sig(j) > 0
      
                delta_Ff(end+1) = double(int( q0+q2*4/L^2*(x*L-x^2), delta_x*(j-1), delta_x*j ));
                delta_Ypsilon(end+1) = delta_Ff(end)/E/A;
                
                
            elseif Sig(j) < 0
        
                delta_Ff(end+1) = -double(int( q0+q2*4/L^2*(x*L-x^2), delta_x*(j-1), delta_x*j ));
                delta_Ypsilon(end+1) = delta_Ff(end)/E/A;
                
                
            else
                delta_Ff(end+1) = 0;
                delta_Ypsilon(end+1) = 0;
            end
      
            
        end
        %calculate the Ypsilon from delta_Ypsilon
        for i = 1:N
            summe=0;
            for j = i:N
                summe = summe+delta_Ypsilon(j);
            end
            Ypsilon(end+1) = summe;
        end
        
        % calculate the u from Ypsilon
        
        for j = 1:N
            summe = 0;
            for s = j:N-1
                summe = summe+(Ypsilon(s)+Ypsilon(s+1))/2*delta_x;
                
            end
            u(end+1) = summe;
        end
        
    end
    u_summe(end+1) = sum(u);
    Ff_summe(end+1) = sum(delta_Ff);
    
end
plot(u_summe(:),Ff_summe(:),'LineWidth',2);
xlabel('u');
ylabel('F');