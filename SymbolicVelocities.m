%----------------------------------------------
%
%       Calculate symbolic Velocities
%
%-----------------------------------------------

clear all, clc

syms q1(t) q2(t) q3(t) q4(t) q5(t) q6(t) q7(t) q8(t) q9(t) a b c d l_femur l_tibia

syms q1dot q2dot q3dot q4dot q5dot q6dot q7dot q8dot q9dot 

PoH = [a;b];
AK = [c;d];

Po = [cos(q1) sin(q1); -sin(q1) cos(q1)]*[-q2; -q3];
H = Po + [cos(q1) sin(q1); -sin(q1) cos(q1)]*PoH;

%right positions
K_r = H + l_femur*[sin(q6-q1); -cos(q6-q1)];
A_r = K_r + l_tibia*[sin(q6-q1+q5); -cos(q6-q1+q5)];

calcn_r = A_r + [cos(q6-q1+q5+q4) -sin(q6-q1+q5+q4); sin(q6-q1+q5+q4) cos(q6-q1+q5+q4)]*AK;
v_calcn_r = diff(calcn_r,t);
v_calcn_r = subs(v_calcn_r,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),diff(q4(t),t),diff(q5(t),t),diff(q6(t),t),diff(q7(t),t),diff(q8(t),t),diff(q9(t),t)],[q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot,]);

toes_r = A_r + [cos(q6-q1+q5+q4) -sin(q6-q1+q5+q4); sin(q6-q1+q5+q4) cos(q6-q1+q5+q4)]*(AK+[0.2;0]);
v_toes_r = diff(toes_r,t);
v_toes_r = subs(v_toes_r,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),diff(q4(t),t),diff(q5(t),t),diff(q6(t),t),diff(q7(t),t),diff(q8(t),t),diff(q9(t),t)],[q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot,]);

%left positions
K_l = H + l_femur*[sin(q7-q1); -cos(q7-q1)];
A_l = K_l + l_tibia*[sin(q7-q1+q8); -cos(q7-q1+q8)];
calcn_l = A_l + [cos(q7-q1+q8+q9) -sin(q7-q1+q8+q9); sin(q7-q1+q8+q9) cos(q7-q1+q8+q9)]*AK;
v_calcn_l = diff(calcn_l,t);
v_calcn_l = subs(v_calcn_l,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),diff(q4(t),t),diff(q5(t),t),diff(q6(t),t),diff(q7(t),t),diff(q8(t),t),diff(q9(t),t)],[q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot,]);

toes_l = A_l + [cos(q7-q1+q8+q9) -sin(q7-q1+q8+q9); sin(q7-q1+q8+q9) cos(q7-q1+q8+q9)]*(AK+[0.2;0]);
v_toes_l = diff(toes_l,t);
v_toes_l = subs(v_toes_l,[diff(q1(t),t),diff(q2(t),t),diff(q3(t),t),diff(q4(t),t),diff(q5(t),t),diff(q6(t),t),diff(q7(t),t),diff(q8(t),t),diff(q9(t),t)],[q1dot, q2dot, q3dot, q4dot, q5dot, q6dot, q7dot, q8dot, q9dot,]);
