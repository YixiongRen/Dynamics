function [Foutx, Fouty, dFdQoutx,dFdQouty, contactflag, contactstatus,totalenergy] = calc_nonlinfricforce_dissi (qamp,kt,fn,mu_x,mu_y,debug,region)
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
dfdrQoldx = 0;   % real part
dfdiQoldx = 0;   % imaginary part
dfdrQoldy = 0;   % real part
dfdiQoldy = 0;   % imaginary part
dfdrQx = zeros(length(qt),1);
dfdiQx = zeros(length(qt),1);
dfdrQy = zeros(length(qt),1);
dfdiQy = zeros(length(qt),1);
statu = [];
contactstatus = 0;

%========================================================================== 
% main iteration loop
%========================================================================== 
for it = 2: length(qt)
    if kt*(qx(it)-x(it-1))>mu_x*fn  % x slip+
        
        
        if kt*(qy(it)-y(it-1))>mu_y*fn %y slip +
            ftx(it) = -mu_x*fn;
            fty(it) = -mu_y*fn;
            dfdrQx(it) = 0;
            dfdiQx(it) = 0;
            dfdrQy(it) = 0;
            dfdiQy(it) = 0;
            x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
            y(it) = qy(it) + fty(it)/kt;
            contactstatus = 11;
            
        elseif kt*(qy(it)-y(it-1))<-mu_y*fn %y slip -
        	ftx(it) = -mu_x*fn;
            fty(it) = +mu_y*fn;
            dfdrQx(it) = 0;
            dfdiQx(it) = 0;
            dfdrQy(it) = 0;
            dfdiQy(it) = 0;
            x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
            y(it) = qy(it) + fty(it)/kt;
            contactstatus = 11;
            
        else %y stick
            ftx(it) = -mu_x*fn;
            x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
            fty(it) = -kt*(qy(it) - y(it-1));
            y(it) = y(it-1);
            dfdrQx(it) = 0;
            dfdiQx(it) = 0;
            dfdrQy(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldy; % the factor 2*2*pi is hard coded due to 2 periods!
            dfdiQy(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldy;
            contactstatus = 10;
            
        end
        
        
        
        
        
    elseif kt*(qx(it)-x(it-1))<-mu_x*fn %x slip -
        
        if kt*(qy(it)-y(it-1))>mu_y*fn %y slip +
            ftx(it) = +mu_x*fn;
            fty(it) = -mu_y*fn;
            dfdrQx(it) = 0;
            dfdiQx(it) = 0;
            dfdrQy(it) = 0;
            dfdiQy(it) = 0;
            x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
            y(it) = qy(it) + fty(it)/kt;
            contactstatus = 11;
            
        elseif kt*(qy(it)-y(it-1))<-mu_y*fn %y slip -
            
            ftx(it) = +mu_x*fn;
            fty(it) = +mu_y*fn;
            dfdrQx(it) = 0;
            dfdiQx(it) = 0;
            dfdrQy(it) = 0;
            dfdiQy(it) = 0;
            x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
            y(it) = qy(it) + fty(it)/kt;
            contactstatus = 11;
        else %y stick
            
            ftx(it) = +mu_x*fn;
            x(it) = qx(it) + ftx(it)/kt; % x == qx- feder-extension
            fty(it) = -kt*(qy(it) - y(it-1));
            y(it) = y(it-1);
            dfdrQx(it) = 0;
            dfdiQx(it) = 0;
            dfdrQy(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldy; % the factor 2*2*pi is hard coded due to 2 periods!
            dfdiQy(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldy;
            contactstatus = 10;
        end        
        
        
        
        
    else  %x stick
        
        if kt*(qy(it)-y(it-1))>mu_y*fn %y slip +
            fty(it) = -mu_y*fn;
            y(it) = qy(it) + fty(it)/kt;
            ftx(it) = -kt*(qx(it) - x(it-1));
            x(it) = x(it-1);
            dfdrQy(it) = 0;
            dfdiQy(it) = 0;
            dfdrQx(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldx; % the factor 2*2*pi is hard coded due to 2 periods!
            dfdiQx(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldx;
            contactstatus = 01;
        elseif kt*(qy(it)-y(it-1))<-mu_y*fn %y slip -
            fty(it) = +mu_y*fn;
            y(it) = qy(it) + fty(it)/kt;
            ftx(it) = -kt*(qx(it) - x(it-1));
            x(it) = x(it-1);
            dfdrQy(it) = 0;
            dfdiQy(it) = 0;
            dfdrQx(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldx; % the factor 2*2*pi is hard coded due to 2 periods!
            dfdiQx(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldx;
            contactstatus = 01;            
            
        else %y stick
            ftx(it) = -kt*(qx(it) - x(it-1));
            fty(it) = -kt*(qy(it) - y(it-1));
            dfdrQx(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldx; % the factor 2*2*pi is hard coded due to 2 periods!
            dfdiQx(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldx;
            dfdrQy(it) = -kt*(cos(2*2*pi/N*(it-1))-cos(2*2*pi/N*(it-2))) + dfdrQoldy; % the factor 2*2*pi is hard coded due to 2 periods!
            dfdiQy(it) = -kt*(sin(2*2*pi/N*(it-1))-sin(2*2*pi/N*(it-2))) + dfdiQoldy;
            x(it) = x(it-1);
            y(it) = y(it-1);
            contactstatus = 00;   
            
            
        end        
        
        
        
    end
    statu(end+1)=contactstatus;
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

Foutx  = [2*real(Fx(length(Fx)/2+2)); -2*imag(Fx(length(Fx)/2+2))]; %first harmonic
Fouty  = [2*real(Fy(length(Fy)/2+2)); -2*imag(Fy(length(Fy)/2+2))]; %first harmonic
dFdQoutx = [2*real(dFdrQx(length(dFdrQx)/2+2)), 2*real(dFdiQx(length(dFdiQx)/2+2)); -2*imag(dFdrQx(length(dFdrQx)/2+2)), -2*imag(dFdiQx(length(dFdiQx)/2+2))];
dFdQouty = [2*real(dFdrQy(length(dFdrQy)/2+2)), 2*real(dFdiQy(length(dFdiQy)/2+2)); -2*imag(dFdrQy(length(dFdrQy)/2+2)), -2*imag(dFdiQy(length(dFdiQy)/2+2))];


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
% Friction Region
if region==1
%     figure(1);
%     subplot(2,2,1);
%     plot(1:1:length(qx),qx,'b--','LineWidth',1);
%     hold on
%     plot(1:1:length(x),x,'r-','LineWidth',1);
%     title('Motion of Contact Pair in x-axis');
%     legend('Disk\_mass','Damper\_mass');
%     xlabel('X Range')
%     ylabel('X axis(m)')
%     grid on
%     subplot(2,2,2);
%     plot(1:1:length(qy),qy,'b--','LineWidth',1);
%     hold on
%     plot(1:1:length(y),y,'r-','LineWidth',1);
%     title('Motion of Contact Pair in y-axis');
%     legend('Disk\_mass','Damper\_mass');
%     xlabel('Y Range')
%     ylabel('Y axis(m)')
%     grid on
%     subplot(2,2,4);
%     
%     leg(1,3) = line([-mu_x*fn/kt,mu_x*fn/kt],[mu_y*fn/kt,mu_y*fn/kt],'LineWidth',2,'LineStyle','--');
%     hold on
%     leg(2,3) = line([-mu_x*fn/kt,-mu_x*fn/kt],[mu_y*fn/kt,-mu_y*fn/kt],'LineWidth',2,'LineStyle','--');
%     hold on
%     leg(3,3) = line([mu_x*fn/kt,mu_x*fn/kt],[mu_y*fn/kt,-mu_y*fn/kt],'LineWidth',2,'LineStyle','--');
%     hold on
%     leg(4,3) = line([-mu_x*fn/kt,mu_x*fn/kt],[-mu_y*fn/kt,-mu_y*fn/kt],'LineWidth',2,'LineStyle','--');
%     hold on
%     statu = [0,statu];
%     for i = length(qx)/2+1:length(qx)
%         if statu(i)==11
%            hold on
%            plot(qx(i)-x(i),qy(i)-y(i),'r.');    
% 
%         elseif statu(i)==10
%            hold on
%             plot(qx(i)-x(i),qy(i)-y(i),'g.');    
%         elseif statu(i)==1
%            hold on
%            plot(qx(i)-x(i),qy(i)-y(i),'b.');             
%         else
%            hold on
%            plot(qx(i)-x(i),qy(i)-y(i),'k.');
% 
%         end
%     end
%     if ishandle(leg(end,1))&& ishandle(leg(end,2))
%         legend([leg(end,1)],'slip',[leg(end,2)],'stick',[leg(end,3)],'sliding boundary');
%     elseif ishandle(leg(end,1))==0 && ishandle(leg(end,2))~=0
%         legend([leg(end,2)],'stick',[leg(end,3)],'sliding boundary');
%     else
%         legend([leg(end,1)],'slip',[leg(end,3)],'sliding boundary');
%     end

%     title('Relative Displacement of Contact Pair in xy-plane');
%     legend('Disk\_mass','Damper\_mass');
%     xlabel('X axis(m)')
%     ylabel('Y axis(m)')
%     grid on
%     axis equal
%     subplot(2,2,3);
%     for i = length(qx)/2+1:length(qx)
%         if statu(i)==11
%            hold on
%            plot(qx(i),qy(i),'r.');
%            hold on
%            plot(x(i),y(i),'rs'); % red for all slip
%         elseif statu(i) == 10
%            hold on
%            plot(qx(i),qy(i),'g.');
%            hold on
%            plot(x(i),y(i),'gs'); % green for x slip y stick
%         elseif statu(i) == 1
%            hold on
%            plot(qx(i),qy(i),'b.');
%            hold on
%            plot(x(i),y(i),'bs'); % blue for x slip y stick
%         else
%            hold on
%            plot(qx(i),qy(i),'k.');
%            hold on
%            plot(x(i),y(i),'ks'); % black for all stick
%         end
%     end
%     title('Motion of Contact Pair in xy-plane');
%     legend('Disk\_mass','Damper\_mass');
%     xlabel('X axis(m)')
%     ylabel('Y axis(m)')
%     grid on
%     axis equal
% %     figure(2)
% %     plot(qx(1:length(x)/2+1),ftx(1:length(ftx)/2+1),'r--');
% %     hold on
% %     plot(qy(1:length(y)/2+1),fty(1:length(fty)/2+1),'b-');
%      
%     figure(666);
%     for i =length(qx)/2+1:20:length(qx)
%         if statu(i)==11
%            hold on
%            plot(qx(i),qy(i),'r.');
%            hold on
%            plot(x(i),y(i),'rs'); % red for all slip
%         elseif statu(i) == 10
%            hold on
%            plot(qx(i),qy(i),'g.');
%            hold on
%            plot(x(i),y(i),'gs'); % green for x slip y stick
%         elseif statu(i) == 1
%            hold on
%            plot(qx(i),qy(i),'b.');
%            hold on
%            plot(x(i),y(i),'bs'); % blue for x slip y stick
%         else
%            hold on
%            plot(qx(i),qy(i),'k.');
%            hold on
%            plot(x(i),y(i),'ks'); % black for all stick
%         end
%     end
%     for i =length(qx)/2+1:100:length(qx)
%         quiver(x(i),y(i),qx(i)-x(i),qy(i)-y(i),'color','r');
%         hold on
%     end
%     title('Motion of Contact Pair in xy-plane');
%     legend('Disk\_mass','Damper\_mass');
%     xlabel('X axis(m)')
%     ylabel('Y axis(m)')
%     grid on
%     axis equal        
%     figure(777);
    dissi_x = [];
    dissi_y = [];
    total_x = [];
    total_y = [];
    for i =length(qx)/2+1:length(qx)-1
        if statu(i)==11
           dissi_x(i)=abs(x(i)-x(i-1))*mu_x*fn;
           dissi_y(i)=abs(y(i)-y(i-1))*mu_y*fn;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;        %  all slip
        elseif statu(i) == 10
           dissi_x(i)=abs(x(i)-x(i-1))*mu_x*fn;
           dissi_y(i)=0;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;      %x slip y stick
        elseif statu(i) == 1
           dissi_x(i)=0;
           dissi_y(i)=abs(y(i)-y(i-1))*mu_y*fn;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;        %  x stick y slip
        else
           dissi_x(i)=0;
           dissi_y(i)=0;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;        %  all stick
        end
    end
    for i=length(qx)/2+1:length(qx)-1
        ratio(i)=(dissi_x(i)+dissi_y(i))/(total_x(i)+total_y(i));
        ratio_x(i)=dissi_x(i)/total_x(i);
        ratio_y(i)=dissi_y(i)/total_y(i);
    end
%     plot(length(qx)/2+1:length(qx)-1,ratio(length(qx)/2+1:length(qx)-1));
%     legend('total energy dissipation ratio[%]');
%     figure(888);
%     plot(length(qx)/2+1:length(qx)-1,ratio_x(length(qx)/2+1:length(qx)-1),'color','b');
%     legend('x direction energy dissipation ratio[%]');
%     figure(999);
%     plot(length(qx)/2+1:length(qx)-1,ratio_y(length(qx)/2+1:length(qx)-1),'color','r');
%     legend('y direction energy dissipation ratio[%]');
%         figure(777);
    dissi_x = [];
    dissi_y = [];
    total_x = [];
    total_y = [];
    for i =2:length(qx)/2
        if statu(i)==11
           dissi_x(i)=abs(x(i)-x(i-1))*mu_x*fn;
           dissi_y(i)=abs(y(i)-y(i-1))*mu_y*fn;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;        %  all slip
        elseif statu(i) == 10
           dissi_x(i)=abs(x(i)-x(i-1))*mu_x*fn;
           dissi_y(i)=0;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;      %x slip y stick
        elseif statu(i) == 1
           dissi_x(i)=0;
           dissi_y(i)=abs(y(i)-y(i-1))*mu_y*fn;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;        %  x stick y slip
        else
           dissi_x(i)=0;
           dissi_y(i)=0;
           total_x(i)=abs(qx(i)-qx(i-1))*mu_x*fn;
           total_y(i)=abs(qy(i)-qy(i-1))*mu_y*fn;        %  all stick
        end
    end
    for i=2:length(qx)/2
        ratio(i)=(dissi_x(i)+dissi_y(i))/(total_x(i)+total_y(i));
        ratio_x(i)=dissi_x(i)/total_x(i);
        ratio_y(i)=dissi_y(i)/total_y(i);
        total(i) = total_x(i)+total_y(i);
        dissi(i) = dissi_x(i)+dissi_y(i);
    end
%     plot(2:length(qx)/2,total(2:length(qx)/2));
%     legend('total energy dissipated [W]');
%     figure(888);
%     plot(2:length(qx)/2,dissi_x(2:length(qx)/2),'color','b');
%     legend('x direction energy dissipated [W]');
%     figure(999);
%     plot(2:length(qx)/2,dissi_y(2:length(qx)/2),'color','r');
%     legend('y direction energy dissipated [W]');
    totalenergy = sum(total(2:length(qx)/2));
else 
    totalenergy = 0;
end
    
end
