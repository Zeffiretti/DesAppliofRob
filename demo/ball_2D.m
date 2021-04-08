function ball

clc
clear all
close all
format long


m = 1.0; 
g = 1;
k_friction=0.2;
k_collision=0.4;

x0=0;
y0=1;
xd0=2.0;
yd0=0;

z0=[x0 y0 xd0 yd0];
t0=0;
dt=10.0;
N_time=10000;
t_ode = t0;
z_ode = z0;
options=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision);

%%%%%%%%%%% bouncing %%%%%%%%%%%%
while 1
    tspan = linspace(t0,t0+dt,N_time);
    [t_temp, z_temp, tfinal] = ode113(@flying_ball,tspan,z0,options,m,g,k_friction,k_collision);
    zplus=collision_ball(t_temp(end),z_temp(end,:),m,g,k_friction,k_collision);
    z0 = zplus;
    t0 = t_temp(end)+dt/N_time;
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    if z_ode(end,4)<power(10,-4)
        break;
    end
end
%%%%%%%%%% slipping %%%%%%%%%%%%%%%%
t0=t_ode(end)+dt/N_time;
z0=z_ode(end,:);
tspan = linspace(t0,t0+dt,N_time);
options=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14);
[t_temp, z_temp] = ode113(@slip_ball,tspan,z0,options,m,g,k_friction,k_collision);
t_ode = [t_ode; t_temp(2:end)];
z_ode = [z_ode; z_temp(2:end,:)];

z=z_ode;
t=t_ode;
% for i=1:4
%     figure(i)
% % plot(t,z(:,i),'ko','LineWidth',1)
% plot(t,z(:,i),'ko')
% end
% figure(5)
% plot(z(:,1),z(:,2),'ko')

fontsize=20;
finalTime = t(end);
animationPlot = plot(z(1,1),z(1,2),'.','MarkerSize',30);
xlabel('x (m)','FontSize',fontsize)
ylabel('y (m)','FontSize',fontsize)

axis([0,40,0,1.5]);
% daspect([1,1,1]);
  
% Animation loop
if 1
currentTime = 0;
tic;
i=1;
while currentTime < finalTime
    
    currentx = interp1(t,z(:,1),currentTime);
    currenty = interp1(t,z(:,2),currentTime);

    plot(currentx,currenty,'ko')
    axis([0,40,0,1.5]);
%     animationPlot.XData = currentx;
%     animationPlot.YData = currenty;
    
    currentTime = toc;
%     drawnow;
    set(gca,'Fontsize',fontsize);
    F(i)=getframe(gcf);
    i=i+1;
end
v = VideoWriter('ball.avi');
open(v);
writeVideo(v,F);
close(v);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zdot = flying_ball(t,z,m,g,k_friction,k_collision)

x=z(1);
y=z(2);
xd=z(3);
yd=z(4);
                                  
xdd=0;
ydd=-g;

zdot = [xd yd xdd ydd]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zplus=collision_ball(t,z,m,g,k_friction,k_collision)      

x=z(1);
y=z(2);
xd=z(3);
yd=z(4);             

zplus = [x y xd -k_collision*yd]; 

function zslip=slip_ball(t,z,m,g,k_friction,k_collision)  
x=z(1);
y=z(2);
xd=z(3);
yd=z(4);
if xd>power(10,-6) || xd<power(10,-6)
    xdd=-xd/abs(xd)*k_friction*g;
else
    xdd=0;
end
zslip = [xd yd xdd 0]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gstop, isterminal,direction]=collision(t,z,m,g,k_friction,k_collision)
x=z(1);
y=z(2);
xd=z(3);
yd=z(4); 
gstop = y;
isterminal=1; %Ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
direction=-1; % The t_final can be approached by any direction is indicated by this