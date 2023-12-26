function PendPendTan = ComputePendPendTan(Xkm1)
% Model equations
import casadi.*

import casadi.*
if isnumeric(Xkm1)
    PendPendTan=zeros(1,12);
else
    PendPendTan=MX.zeros(1,12);
end

q1 = Xkm1(1);
q1dot = Xkm1(2);
q2 = Xkm1(3);
q2dot = Xkm1(4);
q3 = Xkm1(5);
q3dot = Xkm1(6);
q4 = Xkm1(7);
q4dot = Xkm1(8);
q5 = Xkm1(9);
q5dot = Xkm1(10);
q6 = Xkm1(11);
q6dot = Xkm1(12);
q7 = Xkm1(13);
q7dot = Xkm1(14);
q8 = Xkm1(15);
q8dot = Xkm1(16);
q9 = Xkm1(17);
q9dot = Xkm1(18);

% Calculate pos and vel of contact points
PoH = [-0.0707;-0.0661];
a = PoH(1); b = PoH(2);
AK = [-0.04877;-0.04195];
c = AK(1); d = AK(2);
l_femur = 0.41; % m
l_tibia = 0.43; % m

% Computing position
Po = [cos(q1) sin(q1); -sin(q1) cos(q1)]*[-q2; -q3];
H = Po + [cos(q1) sin(q1); -sin(q1) cos(q1)]*PoH;

%right positions
K_r = H +l_femur*[sin(q6-q1); -cos(q6-q1)];
A_r = K_r + l_tibia*[sin(q6-q1+q5); -cos(q6-q1+q5)];
calcn_r = A_r + [cos(q6-q1+q5+q4) -sin(q6-q1+q5-q4); sin(q6-q1+q5-q4) cos(q6-q1+q5-q4)]*AK;
toes_r = A_r + [cos(q6-q1+q5+q4) -sin(q6-q1+q5-q4); sin(q6-q1+q5-q4) cos(q6-q1+q5-q4)]*(AK+[0.2;0]);
v_calcn_r = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - c*sin(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) - q2dot*cos(q1) - l_femur*cos(q1 - q6)*(q1dot - q6dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q5 - q1 + q6)*(q5dot - q1dot + q6dot) - d*cos(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) + q1dot*sin(q1)*q3 + l_femur*sin(q1 - q6)*(q1dot - q6dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q5 - q1 + q6)*(q5dot - q1dot + q6dot) + c*cos(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) + q1dot*cos(q1)*q2];
v_toes_r = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - q2dot*cos(q1) - l_femur*cos(q1 - q6)*(q1dot - q6dot) - sin(q4 - q1 + q5 + q6)*(c + 1/5)*(q4dot - q1dot + q5dot + q6dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q5 - q1 + q6)*(q5dot - q1dot + q6dot) - d*cos(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q4 - q1 + q5 + q6)*(q4dot - q1dot + q5dot + q6dot) + q1dot*sin(q1)*q3 + cos(q4 - q1 + q5 + q6)*(c + 1/5)*(q4dot - q1dot + q5dot + q6dot) + l_femur*sin(q1 - q6)*(q1dot - q6dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q5 - q1 + q6)*(q5dot - q1dot + q6dot) + q1dot*cos(q1)*q2];


%left positions
K_l = H +l_femur*[sin(q7-q1); -cos(q7-q1)];
A_l = K_l + l_tibia*[sin(q7-q1+q8); -cos(q7-q1+q8)];
calcn_l = A_l + [cos(q7-q1+q8+q9) -sin(q7-q1+q8+q9); sin(q7-q1+q8+q9) cos(q7-q1+q8+q9)]*AK;
toes_l = A_l + [cos(q7-q1+q8+q9) -sin(q7-q1+q8+q9); sin(q7-q1+q8+q9) cos(q7-q1+q8+q9)]*(AK+[0.2;0]);
v_calcn_l = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - c*sin(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) - q2dot*cos(q1) - l_femur*cos(q1 - q7)*(q1dot - q7dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q7 - q1 + q8)*(q7dot - q1dot + q8dot) - d*cos(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) + q1dot*sin(q1)*q3 + l_femur*sin(q1 - q7)*(q1dot - q7dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q7 - q1 + q8)*(q7dot - q1dot + q8dot) + c*cos(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) + q1dot*cos(q1)*q2];
v_toes_l = [q1dot*sin(q1)*q2 - q3dot*sin(q1) - q2dot*cos(q1) - l_femur*cos(q1 - q7)*(q1dot - q7dot) - sin(q7 - q1 + q8 + q9)*(c + 1/5)*(q7dot - q1dot + q8dot + q9dot) + b*q1dot*cos(q1) - a*q1dot*sin(q1) + l_tibia*cos(q7 - q1 + q8)*(q7dot - q1dot + q8dot) - d*cos(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) - q1dot*cos(q1)*q3; q2dot*sin(q1) - q3dot*cos(q1) - d*sin(q7 - q1 + q8 + q9)*(q7dot - q1dot + q8dot + q9dot) + q1dot*sin(q1)*q3 + cos(q7 - q1 + q8 + q9)*(c + 1/5)*(q7dot - q1dot + q8dot + q9dot) + l_femur*sin(q1 - q7)*(q1dot - q7dot) - a*q1dot*cos(q1) - b*q1dot*sin(q1) + l_tibia*sin(q7 - q1 + q8)*(q7dot - q1dot + q8dot) + q1dot*cos(q1)*q2];


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

PendPendTan(1) = pen1;
PendPendTan(2) = dpen1;
PendPendTan(3) = dtan1;

PendPendTan(4) = pen2;
PendPendTan(5) = dpen2;
PendPendTan(6) = dtan2;

PendPendTan(7) = pen3;
PendPendTan(8) = dpen3;
PendPendTan(9) = dtan3;

PendPendTan(10) = pen4;
PendPendTan(11) = dpen4;
PendPendTan(12) = dtan4;

end

