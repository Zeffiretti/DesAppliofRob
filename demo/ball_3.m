function ball
clc    %清除原有的记录以及变量
clear all
close all
format long     %设置变量的有效位数为16位
m = 1.0;    %球的质量
g = 9.8;    %重力加速度
k = 1000;    %弹簧刚度
c = 1;    %阻尼系数
k_collision=0.8;    %削弱系数
y0=1;   %垂直初始位置
yd0=0;  %垂直初始速度
z0=[y0 yd0];    %初始状态
t0=0;   %初始时间
dt=10.0;    %仿真时间
N_time=10000;
t_ode = t0;
z_ode = z0;
options_1=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision_1);    %设置ODE参数，设置跳出ODE的events
options_2=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision_2);    %设置ODE参数，设置跳出ODE的events

% 主循环
while 1
    tspan = linspace(t0,t0+dt,N_time);  %t_span代表离散的仿真的节点，从初始到结束
    [t_temp, z_temp, tfinal] = ode113(@flying_ball_1,tspan,z0,options_1,m,g,k_collision,k,c);   %自由落体
    z0 = z_temp(end,:);
    t0 = t_temp(end)+dt/N_time;
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    tspan = linspace(t0,t0+dt,N_time);  %t_span代表离散的仿真的节点，从初始到结束
    [t_temp, z_temp, tfina2] = ode113(@flying_ball_2,tspan,z0,options_2,m,g,k_collision,k,c);   %弹簧阻尼
    z0 = z_temp(end,:);
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
        current_y = interp1(t,z(:,1),currentTime);
        xx=linspace(-1.0,1.0,100);
        yy=linspace(0.0,0.0,100);
        plot(xx,yy,'r-',0,current_y,'ko');
        axis([-1.0,1.0,-0.5,1.5]);
        currentTime = toc;
        set(gca,'Fontsize',fontsize);
        F(i)=getframe(gcf);
        i=i+1;
    end
    % 存储为视频
    v = VideoWriter('ball_k_c.avi');
    open(v);
    writeVideo(v,F);
    close(v);
end


% 子函数——小球在空中对应的状态一阶导
function zdot = flying_ball_1(t,z,m,g,k_collision,k,c)
yd=z(2);                                
ydd=-g;
zdot = [yd ydd]';


% 子函数——小球在弹簧阻尼状态的状态一阶导
function zdot = flying_ball_2(t,z,m,g,k_collision,k,c)
y=z(1);
yd=z(2);                                
ydd=-g-k*y/m-c*yd/m;
zdot=[yd ydd]';


% 子函数——接触弹簧阻尼碰撞设置
function [gstop, isterminal,direction]=collision_1(t,z,m,g,k_collision,k,c)
gstop = z(1);
isterminal=1; 
direction=-1; 


% 子函数——脱离弹簧阻尼碰撞设置
function [gstop, isterminal,direction]=collision_2(t,z,m,g,k_collision,k,c)
gstop = z(1);
isterminal=1; 
direction=1;
