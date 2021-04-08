function ball
clc    %清除原有的记录以及变量
clear all
close all
format long     %设置变量的有效位数为16位
m = 1.0;    %球的质量
g = 9.8;    %重力加速度
k_collision=0.8;    %削弱系数
y0=1;   %垂直初始位置
yd0=0;  %垂直初始速度
x0=0;   %水平初始位置
xd0=0.2;  %水平初始速度
z0=[y0 yd0 x0 xd0];    %初始状态
t0=0;   %初始时间
dt=10.0;    %仿真时间
N_time=10000;
t_ode = t0;
z_ode = z0;
options=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision);    %设置ODE参数，设置跳出ODE的events


% 主循环
while 1
    tspan = linspace(t0,t0+dt,N_time);  %t_span代表离散的仿真的节点，从初始到结束
    [t_temp, z_temp, tfinal] = ode113(@flying_ball,tspan,z0,options,m,g,k_collision);   %自由落体，遇到碰撞跳出
    zplus=collision_ball(t_temp(end),z_temp(end,:),m,g,k_collision);    %碰撞的状态切换
    z0 = zplus;
    t0 = t_temp(end)+dt/N_time;
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    if z_ode(end,2) < power(10,-4)    %一直到碰撞后的速度足够小，结束计算
        break;
    end
end


% 画动态图
z=z_ode;
t=t_ode;
fontsize=20;
finalTime = t(end);
if 1
    currentTime = 0;
    tic;
    i=1;
    while currentTime < finalTime
        current_x = interp1(t,z(:,3),currentTime);
        current_y = interp1(t,z(:,1),currentTime);
        plot(current_x,current_y,'ko');
        axis([-0.5,1.5,0.0,1.5]);
        currentTime = toc;
        set(gca,'Fontsize',fontsize);
        F(i)=getframe(gcf);
        i=i+1;
    end
    % 存储为视频
    v = VideoWriter('ball_two.avi');
    open(v);
    writeVideo(v,F);
    close(v);
end


% 子函数——小球在空中对应的状态一阶导
function zdot = flying_ball(t,z,m,g,k_collision)
yd=z(2);                                
ydd=-g;
xd=z(4);
xdd=0;
zdot = [yd ydd xd xdd]';


% 子函数——碰撞后的状态
function zplus=collision_ball(t,z,m,g,k_collision)      
y=z(1);
yd=z(2);  
x=z(3);
xd=z(4);
zplus = [y -k_collision*yd x xd]; 


% 子函数——碰撞设置
function [gstop, isterminal,direction]=collision(t,z,m,g,k_collision)
y=z(1);
gstop = y;
isterminal=1; %Ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
direction=-1; % The t_final can be approached by any direction is indicated by this