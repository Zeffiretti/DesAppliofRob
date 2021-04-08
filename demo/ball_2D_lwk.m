clc
clear all
close all
format long


m = 1.0; %球的质量
g = 9.8;%重力加速度

k_collision=0.8;
f=0.8;

y0=1;%初始位置
yd0=0;%初始速度
x0=0;
xd0=1;

z0=[y0 yd0 x0 xd0];%初始状态
t0=0;%初始时间
dt=10.0;%仿真时间
N_time=10000;
t_ode = t0;
z_ode = z0;
options=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision);%设置ODE参数，设置跳出ODE的events

%%%%%%%%%%% bouncing %%%%%%%%%%%%
while 1
    tspan = linspace(t0,t0+dt,N_time);
    [t_temp, z_temp, tfinal] = ode113(@flying_ball,tspan,z0,options,m,g,k_collision);%自由落体，遇到碰撞跳出
    zplus=collision_ball(t_temp(end),z_temp(end,:),m,g,k_collision,f);%碰撞的状态切换
    z0 = zplus;
    t0 = t_temp(end)+dt/N_time;
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    if z_ode(end,2)<power(10,-4)%一直到碰撞后的速度足够小，结束计算
        break;
    end
end

z=z_ode;
t=t_ode;
fontsize=20;
nn=length(z(:,1));
%  可视化
    v = VideoWriter('ball.avi');
    open(v);
for i=1:30:nn
    currenty = z(i,1);
    currentx = z(i,3);
    plot(currentx,currenty,'ko');
    axis([-1,5,-0.02,1.5]);
    set(gca,'Fontsize',fontsize);
    F=getframe(gcf);
    
    writeVideo(v,F);

end
    close(v);
energy=m*g*z(:,1)+0.5*m*(z(:,2).^2+z(:,4).^2);

plot(t,energy);
%%%%%%%%%%%%%%微分状态方程%%%%%%%%%%%%%%%%%

function zdot = flying_ball(t,z,m,g,k_collision)
y=z(1);
yd=z(2);                                
ydd=-g;

x=z(3);
xd=z(4);
xdd=0;

zdot = [yd ydd xd xdd]';
end
%%%%%%%%%碰撞的状态转换%%%%%%%%%%

function zplus=collision_ball(t,z,m,g,k_collision,f)      
y=z(1);
yd=z(2); 
x=z(3);
xd=z(4);
zplus = [y -k_collision*yd x f*xd]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gstop, isterminal,direction]=collision(t,z,m,g,k_collision)
y=z(1);
yd=z(2); 
x=z(3);
xd=z(4);
gstop = y;
isterminal=1; %Ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify事件起作用与否
direction=-1; % The t_final can be approached by any direction is indicated by this
end