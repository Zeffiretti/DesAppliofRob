model.L1=0.5;%length of link1
model.L2=0.5;%length of link2//����˼�����

w=pi/2;%��ת���ٶ�
t1=2;%t1֮ǰ��ֱ�ߣ�֮��Ȧ
j=1;
for t=0:0.2:8
    pause(0.01)
    if t<t1
        [x,xd,xdd]=TSpline(0,0,0,0,t1/2,0,0.25*w,t1,t);
        [y,yd,ydd]=TSpline(-1,0,0,-0.9,t1/2,-0.75,0,t1,t);%����������ֵ���滮ĩ�˵��λ�á��ٶȡ����ٶ�
       % z=-1+0.25*t;
       % zd=0.25;
    else
        x=0.25*sin(w*(t-t1));
        xd=0.25*w*cos(w*(t-t1));
        xdd=-0.25*w*w*sin(w*(t-t1));
        y=-0.5-0.25*cos(w*(t-t1));
        yd=0.25*w*sin(w*(t-t1));
        ydd=0.25*w*w*cos(w*(t-t1));%�����Һ����滮ĩ�˵��λ�á��ٶȡ����ٶ�
    end
q2=acos((x*x+y*y-model.L1*model.L1-model.L2*model.L2)/(2*model.L1*model.L2));
q1=atan2(x,-y)-acos((model.L1*model.L1+x*x+y*y-model.L2*model.L2)/(2*model.L1*sqrt(x*x+y*y)));%���˶�ѧ��ؽڽǶ�
% q2=-acos((x*x+y*y-model.L1*model.L1-model.L2*model.L2)/(2*model.L1*model.L2));
% q1=atan2(x,-y)+acos((model.L1*model.L1+x*x+y*y-model.L2*model.L2)/(2*model.L1*sqrt(x*x+y*y)));%���˶�ѧ����һ���
model.Refq(1)=q1;
model.Refq(2)=q2;   
x1=model.L1*sin(q1);
y1=-model.L1*cos(q1);%��L1��ĩ�˵��λ��

plot([0 x1],[0 y1],'r-')
hold on
plot([x1 x],[y1 y],'b-')
hold on
plot(x,y,'ro')
hold on
axis equal;
% plot(y,z,'*')
% hold off
% M(j)=getframe;
% j=j+1;
end
%  movie(M)
