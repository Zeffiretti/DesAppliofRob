function [KE, PE]=energy(t,z,m1, m2, I1, I2, l, a, g)   

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3);I2_1 = I2(1,1);I2_2 = I2(2,2);I2_3 = I2(3,3); 
I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);I2_12 = I2(1,2);I2_13 = I2(1,3);I2_23 = I2(2,3); 
KE = (m2*((a*u2*cos(q1)^2*sin(q2) + a*u2*sin(q1)^2*sin(q2))^2 + (l*u1*cos(q1) - a*u2*cos(q1)*cos(q2) + a*u1*sin(q1)*sin(q2))^2 + (a*u1*cos(q1)*sin(q2) - l*u1*sin(q1) + a*u2*cos(q2)*sin(q1))^2))/2 + u1*((u1*(cos(q2)*(I2_3*cos(q2) + I2_23*sin(q2)) + sin(q2)*(I2_23*cos(q2) + I2_2*sin(q2))))/2 + (u2*sin(q1)*(cos(q2)*(I2_13*sin(q1) + I2_23*cos(q1)*cos(q2) - I2_3*cos(q1)*sin(q2)) + sin(q2)*(I2_12*sin(q1) + I2_2*cos(q1)*cos(q2) - I2_23*cos(q1)*sin(q2))))/2 + (u2*cos(q1)*(cos(q2)*(I2_13*cos(q1) - I2_23*cos(q2)*sin(q1) + I2_3*sin(q1)*sin(q2)) + sin(q2)*(I2_12*cos(q1) - I2_2*cos(q2)*sin(q1) + I2_23*sin(q1)*sin(q2))))/2) + (I1_3*u1^2)/2 + (m1*(a^2*u1^2*cos(q1)^2 + a^2*u1^2*sin(q1)^2))/2 + u2*sin(q1)*((u1*(sin(q1)*(I2_13*cos(q2) + I2_12*sin(q2)) + cos(q1)*cos(q2)*(I2_23*cos(q2) + I2_2*sin(q2)) - cos(q1)*sin(q2)*(I2_3*cos(q2) + I2_23*sin(q2))))/2 + (u2*sin(q1)*(sin(q1)*(I2_1*sin(q1) + I2_12*cos(q1)*cos(q2) - I2_13*cos(q1)*sin(q2)) + cos(q1)*cos(q2)*(I2_12*sin(q1) + I2_2*cos(q1)*cos(q2) - I2_23*cos(q1)*sin(q2)) - cos(q1)*sin(q2)*(I2_13*sin(q1) + I2_23*cos(q1)*cos(q2) - I2_3*cos(q1)*sin(q2))))/2 + (u2*cos(q1)*(sin(q1)*(I2_1*cos(q1) - I2_12*cos(q2)*sin(q1) + I2_13*sin(q1)*sin(q2)) + cos(q1)*cos(q2)*(I2_12*cos(q1) - I2_2*cos(q2)*sin(q1) + I2_23*sin(q1)*sin(q2)) - cos(q1)*sin(q2)*(I2_13*cos(q1) - I2_23*cos(q2)*sin(q1) + I2_3*sin(q1)*sin(q2))))/2) + u2*cos(q1)*((u1*(cos(q1)*(I2_13*cos(q2) + I2_12*sin(q2)) - cos(q2)*sin(q1)*(I2_23*cos(q2) + I2_2*sin(q2)) + sin(q1)*sin(q2)*(I2_3*cos(q2) + I2_23*sin(q2))))/2 + (u2*sin(q1)*(cos(q1)*(I2_1*sin(q1) + I2_12*cos(q1)*cos(q2) - I2_13*cos(q1)*sin(q2)) - cos(q2)*sin(q1)*(I2_12*sin(q1) + I2_2*cos(q1)*cos(q2) - I2_23*cos(q1)*sin(q2)) + sin(q1)*sin(q2)*(I2_13*sin(q1) + I2_23*cos(q1)*cos(q2) - I2_3*cos(q1)*sin(q2))))/2 + (u2*cos(q1)*(cos(q1)*(I2_1*cos(q1) - I2_12*cos(q2)*sin(q1) + I2_13*sin(q1)*sin(q2)) - cos(q2)*sin(q1)*(I2_12*cos(q1) - I2_2*cos(q2)*sin(q1) + I2_23*sin(q1)*sin(q2)) + sin(q1)*sin(q2)*(I2_13*cos(q1) - I2_23*cos(q2)*sin(q1) + I2_3*sin(q1)*sin(q2))))/2);
PE = a*g*m2*cos(q2);