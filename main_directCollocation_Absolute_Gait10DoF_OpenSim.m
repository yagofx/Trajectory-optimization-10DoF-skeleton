function main_directCollocation_Absolute_Gait10DoF_OpenSim

close all,  clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Import CasADi and definition of Implicit or Explicit formulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to where CasADi folder is located
%path('root_to_folder');
import casadi.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Collocation points and definition of constant matrices %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = collocation_points(d, 'legendre');
tau_root=[0 tau_root];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Set constant parameter values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time horizon
T = 2;
namefile = 'motion_L1L2L3-1m2s.mot';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of coordinates and torque controls
nq=10;
ncT=7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Set model variables and model equations  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare model variables
q1 = MX.sym('q1'); % pelvis tilt
q2 = MX.sym('q2'); % pelvis translation x
q3 = MX.sym('q3'); % pelvis translation y
q4 = MX.sym('q4'); % ankle angle right
q5 = MX.sym('q5'); % knee angle right
q6 = MX.sym('q6'); % hip flexion right
q7 = MX.sym('q7'); % hip flexion left
q8 = MX.sym('q8'); % knee angle left
q9 = MX.sym('q9'); % ankle angle left
q10 = MX.sym('q10'); % lumbar extension
q1dot = MX.sym('q1dot');
q2dot = MX.sym('q2dot');
q3dot = MX.sym('q3dot');
q4dot = MX.sym('q4dot');
q5dot = MX.sym('q5dot');
q6dot = MX.sym('q6dot');
q7dot = MX.sym('q7dot');
q8dot = MX.sym('q8dot');
q9dot = MX.sym('q9dot');
q10dot = MX.sym('q10dot');
x=[q1 q1dot q2 q2dot q3 q3dot q4 q4dot q5 q5dot q6 q6dot q7 q7dot q8 q8dot q9 q9dot q10 q10dot];

% uT1 = MX.sym('uT1');
% uT2 = MX.sym('uT2');
% uT3 = MX.sym('uT3');
uT4 = MX.sym('uT4');
uT5 = MX.sym('uT5');
uT6 = MX.sym('uT6');
uT7 = MX.sym('uT7');
uT8 = MX.sym('uT8');
uT9 = MX.sym('uT9');
uT10 = MX.sym('uT10');

uT=[uT4 uT5 uT6 uT7 uT8 uT9 uT10];

ua_q1=MX.sym('ua_q1');
ua_q2=MX.sym('ua_q2');
ua_q3=MX.sym('ua_q3');
ua_q4=MX.sym('ua_q4');
ua_q5=MX.sym('ua_q5');
ua_q6=MX.sym('ua_q6');
ua_q7=MX.sym('ua_q7');
ua_q8=MX.sym('ua_q8');
ua_q9=MX.sym('ua_q9');
ua_q10=MX.sym('ua_q10');
ua=[ua_q1 ua_q2 ua_q3 ua_q4 ua_q5 ua_q6 ua_q7 ua_q8 ua_q9 ua_q10];

% Call Computation Contact Forces
CF = MX.zeros(1,8);

F=external('F','gait10dof.dll');
xdot = [q1dot; ua_q1; q2dot; ua_q2; q3dot; ua_q3; q4dot; ua_q4; q5dot; ua_q5; q6dot; ua_q6; q7dot; ua_q7; q8dot; ua_q8; q9dot; ua_q9; q10dot; ua_q10];

% Cost function
L = sum(uT.^2)+sum(ua.^2)+sum(uT.*(x(8:2:nq*2)));
%+sum(ua.^2)+sum(uT.*(x(8:2:nq*2)));
%sum(x(2:2:nq*2).^2) %angular velocities
%sum(ua.^2); % angular accelerations
%sum(uT.^2) % joint moments o torques
%sum(uT.*(x(8:2:nq*2))); % power consumed

% Boundary conditions
%fixed initial states
initial_states_lb=[0; 0; 0; 0; -1.000; 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0];
initial_states_ub=[0; 0; 0; 0; -0.900; 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0];
initial_states_ig=[0; 0; 0; 0; -0.940; 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0];
%fixed final states
final_states_lb=[0; 0; -1; 0; -1.000; 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0];
final_states_ub=[0; 0; -1; 0; -0.900; 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0];
final_states_ig=[0; 0; -1; 0; -0.940; 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0; 0*(2*pi/360); 0];
%state bounds during the movement

omega_limit = (pi/2)/0.5; % 90º in 0.5 sec
v_limit = 1/0.5; % 1 meter in 0.5 seconds

states_lb=[-5*(2*pi/360); -omega_limit; -2; -v_limit; -2; -v_limit; -pi/2; -omega_limit; -120*(2*pi/360); -omega_limit; -120*(2*pi/360); -omega_limit; -120*(2*pi/360); -omega_limit; -120*(2*pi/360); -omega_limit; -pi/2; -omega_limit; -5*(2*pi/360); -omega_limit]; 
states_ub=[pi/2;   omega_limit;  0.2;  v_limit;  0.2;  v_limit;  pi/2;  omega_limit;  10*(2*pi/360);   omega_limit;  120*(2*pi/360);  omega_limit;  120*(2*pi/360);  omega_limit;  10*(2*pi/360);   omega_limit;  pi/2; omega_limit;   pi/2;  omega_limit];
states_ig=zeros(nq*2,1);
%control bounds
controls_lb=-1000*ones(ncT,1);
controls_ub= 1000*ones(ncT,1);
controls_ig=zeros(ncT,1);
%acceleration bounds (implicit form)
accelerations_lb=-500*ones(nq,1);
accelerations_ub= 500*ones(nq,1);
accelerations_ig= zeros(nq,1);

% Control discretization
N = 20; % number of control intervals (20 intervals/second reasonable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Continuous time dynamics
f = Function('f', {x, uT, ua}, {xdot, L});

h = T/N;

%% time grid
tgrid = linspace(0, T, N+1); % tgrid only contains mesh points
for i=1:4
    dtime(i)=tau_root(i)*T/N;
end

for i=1:N
    tgrid_ext([((i-1)*4+1):1:i*4])=[tgrid(i)+dtime];
end
tgrid_ext(end+1)=T; %tgrid_ext contains mesh points and collocation points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Definition and solve of NLP %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
X0 = MX.sym('X0',nq*2);
w = {w{:}, X0};

lbw = [lbw; initial_states_lb];      %lower bounds of the initial state
ubw = [ubw; initial_states_ub];     %upper bounds of the initial state
w0 =  [ w0; initial_states_ig];      %initial guess of initial state


%% Define collocation points for the states at the first interval, their bounds and initial conditions
Xkm1j={};
for j=1:d
    Xkm1j{j} = MX.sym(['X_0_' num2str(j)], nq*2);
    w = {w{:}, Xkm1j{j}};
    lbw = [lbw; states_lb];         %lower bound for x_kj (j collocation point of x_k-1)    %%%%%%%%%%%%%%%%%%%%%%%%%
    ubw = [ubw; states_ub];         %upper bound for x_kj  (j collocation point of x_k-1)   %%%%%%%%%%%%%%%%%%%%%%%%%
    w0 = [w0; states_ig];   
end

%% Define controls at the first point
Ukm1= MX.sym(['U_0'],ncT);   %torque controls
w = {w{:}, Ukm1};   
lbw = [lbw; controls_lb];       %lower bound for u_k
ubw = [ubw; controls_ub];       %upper bound for u_k
w0 = [w0; controls_ig];

Uakm1 = MX.sym(['Ua_0'],nq);    %acceleration controls
w = {w{:}, Uakm1};    
lbw = [lbw; accelerations_lb];       %lower bound for u_k
ubw = [ubw; accelerations_ub];       %upper bound for u_k
w0 =  [w0;  accelerations_ig];

% Formulate the NLP
Xkm1 = X0;

for k=1:N-1
    % New NLP variable for states
    Xk = MX.sym(['X_' num2str(k)], nq*2);
    w = {w{:}, Xk};
    lbw = [lbw;  states_lb];             
    ubw = [ubw;  states_ub];              
    w0 = [w0; states_ig];
    
    % New NLP variable for the torque control
    Uk = MX.sym(['U_' num2str(k)],ncT);
    w = {w{:}, Uk};   
    lbw = [lbw; controls_lb];        %lower bound for u_k
    ubw = [ubw; controls_ub];        %upper bound for u_k
    w0 = [w0; controls_ig];          %initial guess for u_k

% New NLP variable for the acceleration control
Uak = MX.sym(['U_a' num2str(k)],nq);
w = {w{:},Uak};
lbw = [lbw; accelerations_lb];     %lower bound for ua_k
ubw = [ubw; accelerations_ub];      %upper bound for ua_k
w0 =  [w0;  accelerations_ig];        %initial guess for ua_k

% State and control at collocation points (control is not a design variable at the collocation points)
    Xkj = {};
    Ukj={};
    Uakj={};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nq*2);
        w = {w{:}, Xkj{j}};
        lbw = [lbw; states_lb];         %lower bound for x_kj (j collocation point of x_k)   %%%%%%%%%%%%%%%%%%%%%%%%%
        ubw = [ubw; states_ub];         %upper bound for x_kj  (j collocation point of x_k)   %%%%%%%%%%%%%%%%%%%%%%%%%
        w0 = [w0; states_ig];             %initial guess for x_kj (j collocation point of x_k)  %%%%%%%%%%%%%%%%%%%%%%%%%
       
        %Interpolate controls with lagrange polynomials at the collocation points
        Ukj{j}=MX(ncT,1);
        Uakj{j}=MX(nq,1);
      
        Ukj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Ukm1 Uk]);
        Uakj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Uakm1 Uak]);
    end
        
    % Loop over collocation points
    Xk_end = D(1)*Xkm1;
    for j=1:d
       % Expression for the state derivative at the collocation point
       xp = C(1,j+1)*Xkm1;
       for r=1:d
           xp = xp + C(r+1,j+1)*Xkm1j{r};
       end
      
       % Append collocation equations

        [fj qj]=f(Xkm1j{j},Ukj{j},Uakj{j});
        g = {g{:}, h*fj - xp};               %build defect constraints
        lbg = [lbg; zeros(2*nq,1)];             %lower bounds for defect constraints
        ubg = [ubg; zeros(2*nq,1)];             %upper bounds for defect constraints
      
       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkm1j{j}; ;     %x_k at the last collocation point
  
       % Add contribution to quadrature function
       J = J + B(j+1)*qj*h;
    end    
    
        % path constraints
        % Call everything for CF?
        PendPenDtan = ComputePendPendTan(Xkm1);
        CF(1:2)= ComputationContactForces(PendPenDtan(1),PendPenDtan(2),PendPenDtan(3));
        CF(3:4)= ComputationContactForces(PendPenDtan(4),PendPenDtan(5),PendPenDtan(6));
        CF(5:6)= ComputationContactForces(PendPenDtan(7),PendPenDtan(8),PendPenDtan(9));
        CF(7:8)= ComputationContactForces(PendPenDtan(10),PendPenDtan(11),PendPenDtan(12));

        fIDj=F([Xkm1; Uakm1; CF']);        %ID call
        JM = [0;0;0;Ukm1];
        g = {g{:} fIDj-JM};  %equations of motion as path constraints at mesh points
        lbg = [lbg; zeros(nq,1)];             %lower bounds for path constraints
        ubg = [ubg; zeros(nq,1)];             %upper bounds for path constraints

    
    % Add equality constraints
    g = {g{:}, (Xk_end-Xk)};                        %equality constraint between consecutive mesh points
    lbg = [lbg; zeros(2*nq,1)];                      %lower bound for equality constraints between consecutive mesh points
    ubg = [ubg; zeros(2*nq,1)];                      %upper bound for equality constraints between consecutive mesh points
    
    %save states and controls of the current mesh point
    Ukm1=Uk;
    Uakm1=Uak;
    Xkm1=Xk;
    Xkm1j=Xkj; 
end
  
% path constraints
% Call everything for CF?
PendPenDtan = ComputePendPendTan(Xkm1);
CF(1:2)= ComputationContactForces(PendPenDtan(1),PendPenDtan(2),PendPenDtan(3));
CF(3:4)= ComputationContactForces(PendPenDtan(4),PendPenDtan(5),PendPenDtan(6));
CF(5:6)= ComputationContactForces(PendPenDtan(7),PendPenDtan(8),PendPenDtan(9));
CF(7:8)= ComputationContactForces(PendPenDtan(10),PendPenDtan(11),PendPenDtan(12));

fIDj=F([Xkm1; Uakm1; CF']);        %ID call
JM = [0;0;0;Ukm1];
g = {g{:}, fIDj-JM};  %equations of motion as path constraints at the last minus one point
lbg = [lbg; zeros(nq,1)];             %lower bounds for path constraints
ubg = [ubg; zeros(nq,1)];             %upper bounds for path constraints

%States and controls at the last point
Xk = MX.sym(['X_' num2str(k+1)], nq*2);
w = {w{:}, Xk};
lbw = [lbw; final_states_lb];              %lower bound for last collocation point of x_k %%%%%%%%%%%%%%%%%%%%%%%%%
ubw = [ubw; final_states_ub];              %upper bound for last collocation point of x_k %%%%%%%%%%%%%%%%%%%%%%%%%
w0 = [w0;   final_states_ig];

Uk = MX.sym(['U_' num2str(k+1)],ncT);
w = {w{:}, Uk};  
lbw = [lbw; controls_lb];        %lower bound for u_k
ubw = [ubw; controls_ub];        %upper bound for u_k
w0 = [w0; controls_ig];          %initial guess for u_k

Uak = MX.sym(['U_a' num2str(k+1)],nq);
w = {w{:},Uak};
lbw = [lbw; accelerations_lb];     %lower bound for ua_k
ubw = [ubw; accelerations_ub];      %upper bound for ua_k
w0 = [w0; accelerations_ig];        %initial guess for ua_k

%Interpolate controls with lagrange polynomials at the collocation points of the last mesh interval
for j=1:d
    Ukj{j}=MX(ncT,1);
    Uakj{j}=MX(nq,1);
    Ukj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Ukm1 Uk]);
    Uakj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Uakm1 Uak]);
end

Xk_end = D(1)*Xkm1;
for j=1:d
   % Expression for the state derivative at the collocation point
   xp = C(1,j+1)*Xkm1;
   for r=1:d
       xp = xp + C(r+1,j+1)*Xkm1j{r};
   end

   % Append collocation equations
   % defect constraints
   [fj qj]=f(Xkm1j{j},Ukj{j},Uakj{j});
   g = {g{:}, h*fj-xp};                 %build defect constraints
   lbg = [lbg; zeros(2*nq,1)];             %lower bounds for defect constraints
   ubg = [ubg; zeros(2*nq,1)];             %upper bounds for defect constraints

   % Add contribution to the end state
   Xk_end = Xk_end + D(j+1)*Xkm1j{j};     %x_k at the last collocation point

   % Add contribution to quadrature function
   J = J + B(j+1)*qj*h;
end  

% path constraints
PendPenDtan = ComputePendPendTan(Xkm1);
CF(1:2)= ComputationContactForces(PendPenDtan(1),PendPenDtan(2),PendPenDtan(3));
CF(3:4)= ComputationContactForces(PendPenDtan(4),PendPenDtan(5),PendPenDtan(6));
CF(5:6)= ComputationContactForces(PendPenDtan(7),PendPenDtan(8),PendPenDtan(9));
CF(7:8)= ComputationContactForces(PendPenDtan(10),PendPenDtan(11),PendPenDtan(12));

fIDj=F([Xkm1; Uakm1; CF']);
JM = [0;0;0;Ukm1];
g = {g{:}, fIDj-JM};  %equations of motion as path constraints at the last minus one point
lbg = [lbg; zeros(nq,1)];             %lower bounds for path constraints
ubg = [ubg; zeros(nq,1)];             %upper bounds for path constraints

g = {g{:}, (Xk_end-Xk)};                        %equality constraint between consecutive mesh points
lbg = [lbg; zeros(2*nq,1)];                      %lower bound for equality constraints between consecutive mesh points
ubg = [ubg; zeros(2*nq,1)]; 
    
%% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));

opt.ipopt.max_iter=4000;
opt.ipopt.hessian_approximation = 'limited-memory';
solver = nlpsol('solver', 'ipopt', prob,opt);

%% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Postprocess data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get values for design variables
x_opt=[];
uT_opt=[];
ua_opt=[];
x_opt(1,1:2*nq)=w_opt(1:2*nq);
x_opt_ext=x_opt;
x_opt_ext=[x_opt_ext; reshape(w_opt((2*nq+1):(2*nq+d*2*nq)),2*nq,d)'];
uT_opt(1,1:ncT)=w_opt((2*nq+d*2*nq+1):(2*nq+d*2*nq+ncT));

ua_opt(1,1:nq)=w_opt((2*nq+d*2*nq+2+1):(2*nq+d*2*nq+2+nq));
nvarxint=2*nq+nq+ncT+2*nq*d;
for i=1:N-1
    x_opt=[x_opt; w_opt((2*nq+d*2*nq+ncT+nq+1 +(i-1)*nvarxint):((2*nq+d*2*nq+ncT+nq+2*nq+(i-1)*nvarxint)))'];
    uT_opt=[uT_opt; w_opt((2*nq+d*2*nq+ncT+nq+2*nq+1 +(i-1)*nvarxint):((2*nq+d*2*nq+ncT+nq+2*nq+ncT +(i-1)*nvarxint )))'];
    ua_opt=[ua_opt; w_opt((2*nq+d*2*nq+ncT+nq+2*nq+ncT+1 +(i-1)*nvarxint):((2*nq+d*2*nq+ncT+nq+2*nq+ncT+nq +(i-1)*nvarxint)))'];
    x_opt_ext=[x_opt_ext; x_opt(end,:); reshape(w_opt((2*nq+d*2*nq+ncT+nq+2*nq+ncT+nq+1 +(i-1)*nvarxint):(2*nq+d*2*nq+ncT+nq+2*nq+ncT+nq+2*nq*d +(i-1)*nvarxint)),2*nq,d)'];
end
x_opt=[x_opt; w_opt((   2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+1):((         2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq)))'];
uT_opt=[uT_opt; w_opt(( 2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+1):(     2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+ncT))'];
ua_opt=[ua_opt; w_opt(( 2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+ncT+1):(   2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+ncT+nq))'];
x_opt_ext=[x_opt_ext; x_opt(end,:)];


uT_opt_ext(1:4:(N*(d+1)+1),:)=uT_opt;
ua_opt_ext(1:4:(N*(d+1)+1),:)=ua_opt;
for i=1:N
    for j=1:d
        for qi=1:ncT;
            uT_opt_ext((i-1)*(d+1)+j+1,qi)=LagrangePoly(tau_root(j+1),[0 1],[uT_opt(i,qi) uT_opt(i+1,qi)]);
            ua_opt_ext((i-1)*(d+1)+j+1,qi)=LagrangePoly(tau_root(j+1),[0 1],[ua_opt(i,qi) ua_opt(i+1,qi)]);
        end
    end
end

%Conversion units to degress
% q coordinates
a = x_opt_ext(:,[1 7 9 11 13 15 17 19])*180/pi;
b = [a(:,1) x_opt_ext(:,[3 5]) a(:,2:8)];

% Plotting coordinates
f1 = figure;
subplot(1,2,1);
title('q Coordinates','Interpreter','latex')
yyaxis left
plot(tgrid_ext,b(:,[1 4 5 6 7 8 9 10]),'LineStyle','-');
ylabel('$\theta [deg]$','Interpreter','latex');
yyaxis right
plot(tgrid_ext,b(:,[2 3]),'LineStyle',"-.");
ylabel('$translation [m]$','Interpreter','latex');
xlim([0 T]);
xlabel('time [s]','Interpreter','latex');
legend('q1','q4','q5','q6','q7','q8','q9','q10','q2','q3');
hold all;

% Plotting velocities
subplot(1,2,2);
title('q Velocities','Interpreter','latex')
yyaxis left
plot(tgrid_ext,x_opt_ext(:,[2 8 10 12 14  16 18 20]),'LineStyle','-');
ylabel('$\dot{\theta} [rad/s]$','Interpreter','latex');
yyaxis right
plot(tgrid_ext,x_opt_ext(:,[4 6]),'LineStyle',"-.");
ylabel('$\dot{t  ranslation}[m/s]$','Interpreter','latex');
xlabel('time [s]','Interpreter','latex');
xlim([0 T]);
legend('q1','q4','q5','q6','q7','q8','q9','q10','q2','q3');
hold all;

f2=figure;
title('Moments','Interpreter','latex')
plot(tgrid_ext,uT_opt_ext,'LineStyle','-');
xlim([0 T]);
xlabel('time [s]','Interpreter','latex');
ylabel('$M [Nm]$','Interpreter','latex');
legend('q4','q5','q6','q7','q8','q9','q10');
hold all;

ContactF = zeros(8,81);

for i=1:81

PendPendTan = ComputePendPendTan(x_opt_ext(i,:));
ContactF(1:2,i) = ComputationContactForces(PendPendTan(1),PendPendTan(2),PendPendTan(3));
ContactF(3:4,i) = ComputationContactForces(PendPendTan(4),PendPendTan(5),PendPendTan(6));
ContactF(5:6,i) = ComputationContactForces(PendPendTan(7),PendPendTan(8),PendPendTan(9));
ContactF(7:8,i) = ComputationContactForces(PendPendTan(10),PendPendTan(11),PendPendTan(12));

end

f3=figure;
plot(tgrid_ext,sqrt(ContactF(1,:).^2+ContactF(2,:).^2));
hold on
xlabel('time (s)','Interpreter','latex')
ylabel('Forces (N)','Interpreter','latex')
plot(tgrid_ext,sqrt(ContactF(3,:).^2+ContactF(4,:).^2));
plot(tgrid_ext,sqrt(ContactF(5,:).^2+ContactF(6,:).^2));
plot(tgrid_ext,sqrt(ContactF(7,:).^2+ContactF(8,:).^2));
hold off
legend('$F_{1}$','$F_{2}$','$F_{3}$','$F_{4}$','Interpreter','latex')
title("Reaction forces","Interpreter","latex")

%keyboard;

% Create mot file
q_str.data=[tgrid_ext' b];
q_str.labels={'time','pelvis_tilt','pelvis_tx','pelvis_ty','ankle_angle_r','knee_angle_r','hip_flexion_r','hip_flexion_l','knee_angle_l','ankle_ankle_l','lumbar_extension'};
write_motionFile(q_str,namefile);

%keyboard;

%% Evaluate sparsity of jacobian and hessian
g = vertcat(g{:});
w = vertcat(w{:});
Jac = jacobian(g, w);
figure(4)
spy(sparse(DM.ones(Jac.sparsity())))

keyboard;