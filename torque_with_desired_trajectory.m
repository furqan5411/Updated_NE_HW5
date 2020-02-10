clear all; clc;close all;

syms t I3 I2 I1

g=9.8;
m3 = 60;
m2 = 40;
m1 = 55;
z0 = [0;0;1];

l1=1;
lc1 = 0.5;
l2 =1 ;
lc2 =0.5;
lc3 = 0.5;

I3 = 10;
I2 = 10;
I1 = 10;

A1 = 3;
A2 = 5;
A3 = 7;


q1(t) = A1*sin(t);
q2(t) = A2*cos(2*t);
d3(t) = A3*sin(3*t);

q1_dot(t) = A1*cos(t);
q2_dot(t) = -2*A2*sin(2*t);
d3_dot(t) = 3*A3*cos(3*t);

q1_2dot(t) = diff(q1_dot,t);
q2_2dot(t) = diff(q2_dot,t);
d3_2dot(t) = diff(d3_dot,t);

R0_1 = [cos(q1(t)),0,sin(q1(t));sin(q1(t)),0,cos(q1(t));0,1,0];
R1_2 = [0,-sin(q2(t)),cos(q2(t));0,cos(q2(t)),sin(q2(t));-1,0,0];
R2_3 = [1,0,0;0,1,0;0,0,1];

P_2dot00 = [0;-g;0];

r1_0 = [l1;0;0];
r1_1 = [lc1;0;0];
r2_1 = [l2;0;0];
r2_2 = [lc2;0;0];
r3_2 = [0;d3;0];
r3_3 = [0;lc3;0];

%Forward Recurrsion

w1 = [0;0;q1_dot(t)];; %omega
w1_dot = [0;0;q1_2dot(t)]; %omega dot

w2 = [0;0;q2_dot(t)+q1_dot(t)];
w2_dot = [0;0;q2_2dot(t)+q1_2dot];

w3 = [0;0;q2_dot+q1_dot(t)];
w3_dot = [0;0;q2_2dot(t)+q1_2dot(t)];


P_2dot_1 = R0_1.'*P_2dot00+cross(w1_dot,r1_0)+cross(w1,cross(w1,r1_0));
P_2dot_c1 = P_2dot_1 +cross(w1_dot,r1_1)+cross(w1,cross(w1,r1_1));

P_2dot_2 = R1_2.'*P_2dot_1+cross(w2_dot,r2_1)+cross(w2,cross(w2,r2_1));
P_2dot_c2 = P_2dot_2 +cross(w2_dot,r2_2)+cross(w2,cross(w2,r2_2));

P_2dot_3 = R2_3.'*(P_2dot_2+d3_2dot(t)*z0)+cross(2*d3_dot*w3,(R2_3.')*z0)+cross(w3_dot,r3_2)+cross(w3,cross(w3,r3_2));
P_2dot_c3 = P_2dot_3 +cross(w3_dot,r3_3)+cross(w3,cross(w3,r3_3));

%backward 
f3_3 = m3*P_2dot_c3;
Tau_3(t) = I3*w3_dot+cross(w3,I3*w3)-cross(f3_3,r3_3);

f2_2 = R2_3*f3_3+m2*P_2dot_c2;
Tau_2(t) = R2_3*Tau_3-cross(f2_2,r2_2)-cross(R2_3*f3_3,r2_1-r2_2)+I2*w2_dot+cross(w2,I2*w2);

f1_1 = R1_2*f2_2+m1*P_2dot_c1;
Tau_1(t) = R1_2*Tau_2-cross(f1_1,r1_1)-cross(R1_2*f2_2,r1_0-r1_1)+I1*w1_dot+cross(w1,I1*w1);

a = Tau_1(t);
b = Tau_2(t);
c = Tau_3(t);


% Set latex interpreter for plots text and legends
set(0, 'DefaultTextInterpreter', 'latex') 
set(0, 'DefaultLegendInterpreter', 'latex')

figure;
subplot(3,1,1);
hold on
% fplot(a(1), [0 2*pi],'LineWidth',1.5,'color','red')
% fplot(a(2), [0 2*pi],'LineWidth',1.5,'color','red')
fplot(a(3), [0 2*pi],'LineWidth',1.5,'color','red') %a(3) means along Z
hold off
title('The torque of 1st joint of axis Z');
xlabel('Time $t$ [s]');
ylabel('Torque $\tau$ [$N.m$]');
grid

subplot(3,1,2);
fplot(b(3), [0 2*pi],'LineWidth',1.5,'color','blue')
title('The torque of 2nd joint');
xlabel('Time $t$ [s]');
ylabel('Torque $\tau$ [$N.m$]');
grid

subplot(3,1,3);
fplot(c(3), [0 2*pi],'LineWidth',1.5,'color','green')
title('The torque of 3rd joint');
xlabel('Time $t$ [s]');
ylabel('Torque $\tau$ [$N.m$]');
grid


