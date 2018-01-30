clc;
clear all;
%% the Figure of indenter
% rx = 2;
% ry = 3;
% xx = -4:0.1:4;
% yy = -6:0.1:6;
% [x,y] = meshgrid(xx,yy);
% z = x.^2/4+y.^2/0.6;
% mesh(x,y,z);
%% Definition of the bristles
nx = 21;
ny = 21;
rx = 0.025;       
ry = 0.0025;      
lx = 0.0003;    
ly = 0.0001;      


E1 = 210000000000;
E2 = 300000000;
v1 = 0.269;
v2 = 0.5;
A = 0.0003;
b = 0.04;           % ??????
E_stern = ((1-v1^2)/E1+(1-v2^2)/E2)^-1;
kn = 1.35*E_stern*A/b;
kt = 1.35*2/3*E_stern*A/b;           % ??????
%% Initialization
N = 500;                     % Normal Force:Newton
miu = 0.14;
theta = 0:0.01:5/2*pi;       % array
u_px = 1.8e-6*sin(theta);      % array
v = zeros(nx,ny,4);         % definition of v matrix
for i = 1:length(v)
    for j = 1:length(v)
        x_i = lx*(2*i/nx-1);
        y_i = ly*(2*j/ny-1);
        h = x_i^2/2/rx+y_i^2/2/ry;
        v(i,j,1) = x_i;
        v(i,j,2) = y_i;
        v(i,j,3) = h;
        v(i,j,4) = 1;
    end
end
figure(1);
mesh(v(:,:,1),v(:,:,2),v(:,:,3));   % to see if v is created correctly

%% Hesteresis
u_pz =1e-06;           % ??????
fn = zeros(nx,ny);
Gap = zeros(nx,ny);
w = v(:,:,1);
ft = zeros(nx,ny);
UPX = [];
FPX = [];
UPX(1) = 0;
FPX(1) = 0;
delta_t = 1e-10;
for l = 2:length(u_px)

    UPX(end+1) = u_px(l); 
    delta_upx = u_px(l) - u_px(l-1);
    velocity = delta_upx/delta_t;  
    miu = 1.4 + 0.5*exp(-9*abs(velocity));
    for i = 1:length(fn)
        for j = 1:length(fn)
            w(i,j) = w(i,j) + delta_upx;            
            Gap(i,j) = v(i,j,3);
            Gap(i,j) = Gap(i,j) - u_pz;
            if Gap(i,j)<0
                fn(i,j) = -kn*Gap(i,j);
                ft(i,j) = kt*(w(i,j) - v(i,j,1));
                if abs(ft(i,j))>=miu*abs(fn(i,j))
                        ft(i,j) = miu*ft(i,j)/abs(ft(i,j))*abs(fn(i,j));    % interesting!!!
                        w(i,j) = v(i,j,1) + ft(i,j)/kt;
                end           
                    
            else
                w(i,j) = v(i,j,1);
            end
        
            
        end
    end
    figure(2);
    quiver(0:4:80,0:1:20,ft,zeros(size(ft)));
    set (gcf,'Position',[600,300,800,400], 'color','w');    
    pause(0.01);
    FPX(end+1) = sum(ft(:));
end

%%
FPX = FPX/100;
figure(3);
h(1)=plot(UPX(1:end/5+1),FPX(1:end/5+1),'r--','LineWidth',1);
hold on
h(2)=plot(UPX(end/5+1:end/5*3+1),FPX(end/5+1:end/5*3+1),'g-','LineWidth',1);
hold on
h(3)=plot(UPX(end/5*3+1:end),FPX(end/5*3+1:end),'b:','LineWidth',1);
legend([h(1) h(2) h(3)],{'$$ 0\leq \theta \le \frac{\pi}{2} $$',...
    '$$ \frac{\pi}{2}\leq \theta \le \frac{3\pi}{2} $$',...
    '$$ \frac{3\pi}{2}\leq \theta \le \frac{5\pi}{2} $$'},'Interpreter','latex','Location','northwest','FontSize',16);
xlabel('u [m]');
ylabel('F [N]');
title('The force response when displacement ux is given as: $$1.8e^{-6}\times\sin(\theta)$$','Interpreter','latex','FontSize',16);
grid on