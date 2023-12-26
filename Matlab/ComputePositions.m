%--------------------------------------------
%
%       Compute penetration, dpen & dtan
%
%--------------------------------------------

import casadi.*

close all,  clc, clear all
% Calculate pos and vel of contact points
PoH = [-0.0707; -0.0661]; % a; b
a = PoH(1); b = PoH(2);
AK = [-0.04877;-0.04195]; % c; d
c = AK(1); d = AK(2);
torso_offest = [-0.1007; 0.0815];

l_femur = 0.41; % m
l_tibia = 0.43; % m

contactSphereRadius=0.01;

% Model variables
q1 = 15.882*(2*pi/360); % pelvis tilt
q1dot = 0;
q2 = -0.5; % pelvis translation x
q2dot = 0;
q3 = -0.976; % pelvis translation y
q3dot = 0;
q4 = -15.88*(2*pi/360); % ankle angle right
q4dot = 0;
q5 = -20*(2*pi/360); % knee angle right
q5dot = -3;
q6 = 26.218*(2*pi/360); % hip flexion right
q6dot = 0;
q7 = -21.176*(2*pi/360); % hip flexion left
q7dot = 0;
q8 = -32.059*(2*pi/360); % knee angle left
q8dot = 0;
q9 = -26.471*(2*pi/360); % ankle angle left
q9dot = 0;


q10 = 0; % lumbar extension

% Computing position
Po = [cos(q1) sin(q1); -sin(q1) cos(q1)]*[-q2; -q3];
H = Po + [cos(q1) sin(q1); -sin(q1) cos(q1)]*PoH;
T = Po + [cos(q1) sin(q1); -sin(q1) cos(q1)]*torso_offest;

%right positions
K_r = H + l_femur*[sin(q6-q1); -cos(q6-q1)];
A_r = K_r + l_tibia*[sin(q6-q1+q5); -cos(q6-q1+q5)];

calcn_r = A_r + [cos(q6-q1+q5+q4) -sin(q6-q1+q5+q4); sin(q6-q1+q5+q4) cos(q6-q1+q5+q4)]*AK;
v_calcn_r = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - c*sin(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) - q2dot*cos(q1) - l_femur*cos(q1 - q6)*(q1dot - q6dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q5 - q1 + q6)*(q5dot - q1dot + q6dot) - d*cos(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) + q1dot*sin(q1)*q3 + l_femur*sin(q1 - q6)*(q1dot - q6dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q5 - q1 + q6)*(q5dot - q1dot + q6dot) + c*cos(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) + q1dot*cos(q1)*q2];

toes_r = A_r + [cos(q6-q1+q5+q4) -sin(q6-q1+q5+q4); sin(q6-q1+q5+q4) cos(q6-q1+q5+q4)]*(AK+[0.2;0]);
v_toes_r = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - q2dot*cos(q1) - l_femur*cos(q1 - q6)*(q1dot - q6dot) - sin(q4 - q1 + q5 + q6)*(c + 1/5)*(q4dot - q1dot + q5dot + q6dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q5 - q1 + q6)*(q5dot - q1dot + q6dot) - d*cos(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) + q1dot*sin(q1)*q3 + cos(q4 - q1 + q5 + q6)*(c + 1/5)*(q4dot - q1dot + q5dot + q6dot) + l_femur*sin(q1 - q6)*(q1dot - q6dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q5 - q1 + q6)*(q5dot - q1dot + q6dot) + q1dot*cos(q1)*q2];
%left positions
K_l = H + l_femur*[sin(q7-q1); -cos(q7-q1)];
A_l = K_l + l_tibia*[sin(q7-q1+q8); -cos(q7-q1+q8)];

calcn_l = A_l + [cos(q7-q1+q8+q9) -sin(q7-q1+q8+q9); sin(q7-q1+q8+q9) cos(q7-q1+q8+q9)]*AK;
v_calcn_l = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - c*sin(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) - q2dot*cos(q1) - l_femur*cos(q1 - q7)*(q1dot - q7dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q7 - q1 + q8)*(q7dot - q1dot + q8dot) - d*cos(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) + q1dot*sin(q1)*q3 + l_femur*sin(q1 - q7)*(q1dot - q7dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q7 - q1 + q8)*(q7dot - q1dot + q8dot) + c*cos(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) + q1dot*cos(q1)*q2];

toes_l = A_l + [cos(q7-q1+q8+q9) -sin(q7-q1+q8+q9); sin(q7-q1+q8+q9) cos(q7-q1+q8+q9)]*(AK+[0.2;0]);
v_toes_l = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - q2dot*cos(q1) - l_femur*cos(q1 - q7)*(q1dot - q7dot) - sin(q7 - q1 + q8 + q9)*(c + 1/5)*(q7dot - q1dot + q8dot + q9dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q7 - q1 + q8)*(q7dot - q1dot + q8dot) - d*cos(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) + q1dot*sin(q1)*q3 + cos(q7 - q1 + q8 + q9)*(c + 1/5)*(q7dot - q1dot + q8dot + q9dot) + l_femur*sin(q1 - q7)*(q1dot - q7dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q7 - q1 + q8)*(q7dot - q1dot + q8dot) + q1dot*cos(q1)*q2];

%Plot the sequence of joint positions of human body
f1 = figure;
plot([T(1) Po(1) H(1) K_r(1) A_r(1) calcn_r(1) toes_r(1)],[T(2) Po(2) H(2) K_r(2) A_r(2) calcn_r(2) toes_r(2)],'LineWidth',2);
circle(calcn_r(1), calcn_r(2),contactSphereRadius);
circle(toes_r(1), toes_r(2),contactSphereRadius);
hold on
plot([T(1) Po(1) H(1) K_l(1) A_l(1) calcn_l(1) toes_l(1)],[T(2) Po(2) H(2) K_l(2) A_l(2) calcn_l(2) toes_l(2)],'LineWidth',2);
circle(calcn_l(1), calcn_l(2),contactSphereRadius);
circle(toes_l(1), toes_l(2),contactSphereRadius);
hold on
plot(0,0,'ok','MarkerSize',8);
plot([1 0 0],[0 0 1],'k');
xlim([-0.5 1.5]); 
ylim([-0.5 1.5]);
hold off

% F1
pen1 = -toes_r(2);
dpen1 = -v_toes_r(2);
dtan1 = v_toes_r(1);
% F2
pen2 = -calcn_r(2);
dpen2 = -v_calcn_r(2);
dtan2 = v_calcn_r(1);
%F3
pen3 = -calcn_l(2);
dpen3 = -v_calcn_l(2);
dtan3 = v_calcn_l(1);
%F4
pen4 = -toes_l(2);
dpen4 = -v_toes_l(2);
dtan4 = v_toes_l(1);

% Call Computation Contact Forces
CF = MX.zeros(1,8);
CF(1:2)= ComputationContactForces(pen1,dpen1,dtan1);
CF(3:4)= ComputationContactForces(pen2,dpen2,dtan2);
CF(5:6)= ComputationContactForces(pen3,dpen3,dtan3);
CF(7:8)= ComputationContactForces(pen4,dpen4,dtan4);


function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
end



