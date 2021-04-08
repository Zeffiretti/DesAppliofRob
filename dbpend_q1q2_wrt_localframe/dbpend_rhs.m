function zdot=dbpend_rhs(t,z,flag,m1, m2, I1, ...
                              I2, l, a, g);           

q1 = z(1);                          
u1 = z(2);                          
q2 = z(3);                         
u2 = z(4);                       
J =[ l*cos(q1 + q2) + l*cos(q1), l*cos(q1 + q2);
      l*sin(q1 + q2) + l*sin(q1), l*sin(q1 + q2)];
t1 = 2;
w=pi/2;
if t<t1
        [x,xd,xdd]=TSpline(0,0,0,0,t1/2,0,0.25*w,t1,t);
        [y,yd,ydd]=TSpline(-1,0,0,-0.9,t1/2,-0.75,0,t1,t);%三次样条差值，规划末端点的位置、速度、加速度
       % z=-1+0.25*t;
       % zd=0.25;
else
        x=0.25*sin(w*(t-t1));
        xd=0.25*w*cos(w*(t-t1));
        xdd=-0.25*w*w*sin(w*(t-t1));
        y=-0.5-0.25*cos(w*(t-t1));
        yd=0.25*w*sin(w*(t-t1));
        ydd=0.25*w*w*cos(w*(t-t1));%用正弦函数规划末端点的位置、速度、加速度
end
q2_ref=acos((x*x+y*y-l*l-l*l)/(2*l*l));
q1_ref=atan2(x,-y)-acos((l*l+x*x+y*y-l*l)/(2*l*sqrt(x*x+y*y)));%逆运动学求关节角度
temp = [xd yd];
u = temp/J;
u1_ref = u(1);
u2_ref = u(2);

% w=2;
% q1_ref = 1-cos(w*t);
% u1_ref = w*sin(w*t);
% q2_ref = 1-cos(w*t);
% u2_ref = w*sin(w*t);
kp = 1000;
kd=100;
T1 = kp*(q1_ref-q1)+kd*(u1_ref-u1); %zero torques for unactuated system
T2 = kp*(q2_ref-q2)+kd*(u2_ref-u2); 
% T1=0;
% T2=0;
M11 = - I1 - I2 - a^2*m1 - a^2*m2 - l^2*m2 - 2*a*l*m2*cos(q2);%-I1-a^2*m1-m2*l^2;
M12 = - m2*a^2 - l*m2*cos(q2)*a - I2;%-cos(-q1+q2)*l*a*m2;
M21 = - m2*a^2 - l*m2*cos(q2)*a - I2;%-cos(-q1+q2)*l*a*m2;
M22 = - m2*a^2 - I2;%-m2*a^2-I2;

RHS1 = - a*l*m2*sin(q2)*u2^2 - 2*a*l*m2*u1*sin(q2)*u2 - T1 + a*g*m2*sin(q1 + q2) + a*g*m1*sin(q1) + g*l*m2*sin(q1);%m2*g*l*sin(q1)-m2*l*u2^2*a*sin(-q1+q2)+a*sin(q1)*m1*g-T1+T2;
RHS2 = a*l*m2*sin(q2)*u1^2 - T2 + a*g*m2*sin(q1 + q2);%a*sin(q2)*m2*g-T2+m2*a*u1^2*l*sin(-q1+q2);

M    = [M11 M12; M21 M22];
RHS  = [RHS1 ; RHS2];
udot =  M \ RHS;

ud1 = udot(1);
ud2 = udot(2);

zdot = [u1 ud1 u2 ud2]'  ;            
