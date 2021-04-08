function ball
clc    %���ԭ�еļ�¼�Լ�����
clear all
close all
format long     %���ñ�������Чλ��Ϊ16λ
m = 1.0;    %�������
g = 9.8;    %�������ٶ�
k_collision=0.8;    %����ϵ��
y0=1;   %��ֱ��ʼλ��
yd0=0;  %��ֱ��ʼ�ٶ�
x0=0;   %ˮƽ��ʼλ��
xd0=0.2;  %ˮƽ��ʼ�ٶ�
z0=[y0 yd0 x0 xd0];    %��ʼ״̬
t0=0;   %��ʼʱ��
dt=10.0;    %����ʱ��
N_time=10000;
t_ode = t0;
z_ode = z0;
options=odeset('abstol',2.25*1e-14,'reltol',2.25*1e-14,'events',@collision);    %����ODE��������������ODE��events


% ��ѭ��
while 1
    tspan = linspace(t0,t0+dt,N_time);  %t_span������ɢ�ķ���Ľڵ㣬�ӳ�ʼ������
    [t_temp, z_temp, tfinal] = ode113(@flying_ball,tspan,z0,options,m,g,k_collision);   %�������壬������ײ����
    zplus=collision_ball(t_temp(end),z_temp(end,:),m,g,k_collision);    %��ײ��״̬�л�
    z0 = zplus;
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
        current_x = interp1(t,z(:,3),currentTime);
        current_y = interp1(t,z(:,1),currentTime);
        plot(current_x,current_y,'ko');
        axis([-0.5,1.5,0.0,1.5]);
        currentTime = toc;
        set(gca,'Fontsize',fontsize);
        F(i)=getframe(gcf);
        i=i+1;
    end
    % �洢Ϊ��Ƶ
    v = VideoWriter('ball_two.avi');
    open(v);
    writeVideo(v,F);
    close(v);
end


% �Ӻ�������С���ڿ��ж�Ӧ��״̬һ�׵�
function zdot = flying_ball(t,z,m,g,k_collision)
yd=z(2);                                
ydd=-g;
xd=z(4);
xdd=0;
zdot = [yd ydd xd xdd]';


% �Ӻ���������ײ���״̬
function zplus=collision_ball(t,z,m,g,k_collision)      
y=z(1);
yd=z(2);  
x=z(3);
xd=z(4);
zplus = [y -k_collision*yd x xd]; 


% �Ӻ���������ײ����
function [gstop, isterminal,direction]=collision(t,z,m,g,k_collision)
y=z(1);
gstop = y;
isterminal=1; %Ode should terminate is conveyed by 1, if you put 0 it goes till the final time u specify
direction=-1; % The t_final can be approached by any direction is indicated by this