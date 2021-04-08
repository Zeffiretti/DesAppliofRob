function ball
clc    %���ԭ�еļ�¼�Լ�����
clear all
close all
format long     %���ñ�������Чλ��Ϊ16λ
m = 1.0;    %�������
g = 9.8;    %�������ٶ�
k = 1000;    %���ɸն�
c = 1;    %����ϵ��
k_collision=0.8;    %����ϵ��
y0=1;   %��ֱ��ʼλ��
yd0=0;  %��ֱ��ʼ�ٶ�
z0=[y0 yd0];    %��ʼ״̬
t0=0;   %��ʼʱ��
dt=10.0;    %����ʱ��
N_time=10000;
t_ode = t0;
z_ode = z0;
options_1=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision_1);    %����ODE��������������ODE��events
options_2=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision_2);    %����ODE��������������ODE��events

% ��ѭ��
while 1
    tspan = linspace(t0,t0+dt,N_time);  %t_span������ɢ�ķ���Ľڵ㣬�ӳ�ʼ������
    [t_temp, z_temp, tfinal] = ode113(@flying_ball_1,tspan,z0,options_1,m,g,k_collision,k,c);   %��������
    z0 = z_temp(end,:);
    t0 = t_temp(end)+dt/N_time;
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    tspan = linspace(t0,t0+dt,N_time);  %t_span������ɢ�ķ���Ľڵ㣬�ӳ�ʼ������
    [t_temp, z_temp, tfina2] = ode113(@flying_ball_2,tspan,z0,options_2,m,g,k_collision,k,c);   %��������
    z0 = z_temp(end,:);
    t0 = t_temp(end)+dt/N_time;
    t_ode = [t_ode; t_temp(2:end); t0];
    z_ode = [z_ode; z_temp(2:end,:); z0];
    if z_ode(end,2) < power(10,-4)    %һֱ����ײ����ٶ��㹻С����������
        break;
    end
end


% ����̬ͼ
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
    % �洢Ϊ��Ƶ
    v = VideoWriter('ball_k_c.avi');
    open(v);
    writeVideo(v,F);
    close(v);
end


% �Ӻ�������С���ڿ��ж�Ӧ��״̬һ�׵�
function zdot = flying_ball_1(t,z,m,g,k_collision,k,c)
yd=z(2);                                
ydd=-g;
zdot = [yd ydd]';


% �Ӻ�������С���ڵ�������״̬��״̬һ�׵�
function zdot = flying_ball_2(t,z,m,g,k_collision,k,c)
y=z(1);
yd=z(2);                                
ydd=-g-k*y/m-c*yd/m;
zdot=[yd ydd]';


% �Ӻ��������Ӵ�����������ײ����
function [gstop, isterminal,direction]=collision_1(t,z,m,g,k_collision,k,c)
gstop = z(1);
isterminal=1; 
direction=-1; 


% �Ӻ����������뵯��������ײ����
function [gstop, isterminal,direction]=collision_2(t,z,m,g,k_collision,k,c)
gstop = z(1);
isterminal=1; 
direction=1;
