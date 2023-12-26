%--------------------------------------------
%
%       Compute CF
%
%--------------------------------------------

close all

ContactF = zeros(8,81);

for i=1:81

PendPendTan = ComputePendPendTan(x_opt_ext(i,:));
ContactF(1:2,i) = ComputationContactForces(PendPendTan(1),PendPendTan(2),PendPendTan(3));
ContactF(3:4,i) = ComputationContactForces(PendPendTan(4),PendPendTan(5),PendPendTan(6));
ContactF(5:6,i) = ComputationContactForces(PendPendTan(7),PendPendTan(8),PendPendTan(9));
ContactF(7:8,i) = ComputationContactForces(PendPendTan(10),PendPendTan(11),PendPendTan(12));

end

figure;
plot(time,sqrt(ContactF(1,:).^2)+(ContactF(2,:).^2));
hold on
xlabel('time (s)','Interpreter','latex')
ylabel('Forces (N)','Interpreter','latex')
plot(time,sqrt(ContactF(3,:).^2)+(ContactF(4,:).^2));
plot(time,sqrt(ContactF(5,:).^2)+(ContactF(6,:).^2));
plot(time,sqrt(ContactF(7,:).^2)+(ContactF(8,:).^2));
hold off
legend('$F_{1}$','$F_{2}$','$F_{3}$','$F_{4}$','Interpreter','latex')
title("Reaction forces","Interpreter","latex")



