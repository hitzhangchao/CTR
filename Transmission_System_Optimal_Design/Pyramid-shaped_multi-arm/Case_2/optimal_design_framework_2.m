%% Optimal Design Framework for Pyramid-shaped Transmission System
% This case demonstrates the design optimization of a four-arm CTR with 3 tubes in each arm.
% 1. Establish the design constraints
% 2. Formulate optimization problem
% 3. Getting the Pareto front
% by Chao Zhang
% Date：2023/10/1

clear;clc;
syms ODs_1 ODs_2 ODs_3 lb alpha le beta_1 beta_2 beta_3 real    % Design variables


%% Load CTR params
CTR_params_2;


%% Geometric Constraints - Relations of tube diameters
IDs_3 = ODn_3+2*b;                      % Inner diameter of SS tube 2 (inner tube)
%ODs_1 <= min(a,D-a/sin(pi/n));
%ODs_1 >= ODn_1+2*b+2*w;
%ODs_2 >= max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w)
%ODs_2 <= ODs_1-2*w-2*delta
%ODs_3 >= max(ODn_2+2*b-2*delta,ODn_3+2*b+2*w)
%ODs_3 <= ODs_2-2*w-2*delta


%% Geometric Constraints - Shape and dimensions of SS tubes
v = 2/3 * lb * alpha;                   % Maximum deflection of SS tube
L1 = delta_L_1+le+lb+ls;                % Length of SS tube 1
L2 = delta_L_2+le+lb+ls-la;             % Length of SS tube 2
L3 = delta_L_3+le+lb+ls-2*la;           % Length of SS tube 3
% le >= le_min


%% Geometric Constraints - Collision avoidance of actuation units
% Coordinates of tube points - Here we directly let beta=90° to accelerate the computing time
O = [0;0;0];
O1 = [sqrt(2)/2*a;0;0];
O2 = [0;sqrt(2)/2*a;0];
O3 = [-sqrt(2)/2*a;0;0];
O4 = [0;-sqrt(2)/2*a;0];

A = [0;0;ls];
A1 = [sqrt(2)/2*a;0;ls];
A2 = [0;sqrt(2)/2*a;ls];
A3 = [-sqrt(2)/2*a;0;ls];
A4 = [0;-sqrt(2)/2*a;ls];

B = [0;0;ls+lb];
B1 = [sqrt(2)/2*a+v;0;ls+lb];
B2 = [0;sqrt(2)/2*a+v;ls+lb];
B3 = [-sqrt(2)/2*a-v;0;ls+lb];
B4 = [0;-sqrt(2)/2*a-v;ls+lb];

C = [0;0;ls+lb+le*cos(alpha)];
C1 = [sqrt(2)/2*a+v+le*sin(alpha);0;ls+lb+le*cos(alpha)];
C2 = [0;sqrt(2)/2*a+v+le*sin(alpha);ls+lb+le*cos(alpha)];
C3 = [-sqrt(2)/2*a-v-le*sin(alpha);0;ls+lb+le*cos(alpha)];
C4 = [0;-sqrt(2)/2*a-v-le*sin(alpha);ls+lb+le*cos(alpha)];

% line vector
C1C2 = C2-C1;
C2C3 = C3-C2;
C3C4 = C4-C3;
C4C1 = C1-C4;
d_C1C2 = sqrt(C1C2(1)^2+C1C2(2)^2);
d_C2C3 = sqrt(C2C3(1)^2+C2C3(2)^2);
d_C3C4 = sqrt(C3C4(1)^2+C3C4(2)^2);
d_C4C1 = sqrt(C4C1(1)^2+C4C1(2)^2);


% Normal vector of the plane \kappa_j
n1 = (C1-B1)/le;
n2 = (C2-B2)/le;
n3 = (C3-B3)/le;
n4 = (C4-B4)/le;

% Included angle between the line and the plane
sin_gama12 = dot(n1,C1C2)/d_C1C2;
sin_gama21 = dot(n2,-C1C2)/d_C1C2;
sin_gama23 = dot(n2,C2C3)/d_C2C3;
sin_gama32 = dot(n3,-C2C3)/d_C2C3;
sin_gama34 = dot(n3,C3C4)/d_C3C4;
sin_gama43 = dot(n4,-C3C4)/d_C3C4;
sin_gama41 = dot(n4,C4C1)/d_C4C1;
sin_gama14 = dot(n1,-C4C1)/d_C4C1;

% Projecting the radius ra
d_C1D12 = ra*sqrt(1-sin_gama12^2);
d_C2D21 = ra*sqrt(1-sin_gama21^2);
d_C2D23 = ra*sqrt(1-sin_gama23^2);
d_C3D32 = ra*sqrt(1-sin_gama32^2);
d_C3D34 = ra*sqrt(1-sin_gama34^2);
d_C4D43 = ra*sqrt(1-sin_gama43^2);
d_C4D41 = ra*sqrt(1-sin_gama41^2);
d_C1D14 = ra*sqrt(1-sin_gama14^2);

% Minimum distance between adjacent actuation units
d1 = d_C1C2 - d_C1D12 - d_C2D21;
d2 = d_C2C3 - d_C2D23 - d_C3D32;
d3 = d_C3C4 - d_C3D34 - d_C4D43;
d4 = d_C4C1 - d_C4D41 - d_C1D14;


%% Deformation Constraints - Elastic strain limit of bending
epsilon = ODs_1*alpha/lb;               % Maximum bendign strain in SS tube
%epsilon <= epsilon_e


%% Deformation Constraints - Torsional twisting of SS tubes
% Inner SS tube
J3 = pi*(ODs_3^4-IDs_3^4)/32;           % Cross-sectional polar moment of inertia, m^4
kz_3 = Gs*J3;                           % Torsion stiffness, m^4*Pa
zeta_3 = tau_3_max*(L3)/kz_3;           % Maximum relative twist angle, rad

% Middle SS tube
IDs_2 = ODs_3+2*delta;
J2 = pi*(ODs_2^4-IDs_2^4)/32;
kz_2 = Gs*J2;
zeta_2 = tau_2_max*(L2)/kz_2;

% Outer SS tube
IDs_1 = ODs_2+2*delta;
J1 = pi*(ODs_1^4-IDs_1^4)/32;
kz_1 = Gs*J1;
zeta_1 = tau_1_max*(L1)/kz_1;


%% function and constraints generation
% !! Update "multiobjective_func.m" and "nonlinear_constraints.m" manually

% F(1) = zeta_i
Func_1 = (180/pi)*zeta_3

% F(2) = le+lb
Func_2 = le + lb + ls - 2*la + delta_L_3

% ODs_3 <= ODs_2-2*w-2*delta
Cons_1 = ODs_3 - ODs_2 + 2*w + 2*delta

% ODs_2 <= ODs_1-2*w-2*delta
Cons_2 = ODs_2 - ODs_1 + 2*w + 2*delta

% epsilon <= epsilon_e
Cons_3 = epsilon - epsilon_e

% zeta_i<=zeta_max
Cons_4 = (180/pi)*(zeta_3 - zeta_max)

% d1>=0 -> -d1<=0 (make beta_1 = pi/2)
Cons_5 = - subs(d1,beta_1,2*pi/n)

% zeta_3=zeta_2
Cons_6 = (180/pi)*(zeta_3 - zeta_2)

% zeta_3=zeta_1
Cons_7 = (180/pi)*(zeta_3 - zeta_1)


%% Formulation of the multiobjective optimization problem
tic
nvars = 6;                              %numbers of optimization variables
func = @multiobjective_func_2;            %objective function
Aineq = [-1,1,0,0,0,0;0,-1,1,0,0,0];                   %linear constraint: Ax<=b
bineq = [-2*w-2*delta;-2*w-2*delta];                   %linear inequally constraint: Ax<=b
Aeq = [];                               %linear equally constraints:Aeq x = beq
beq = [];
%lbd = [ODn_1+2*b+2*w,max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),max(ODn_2+2*b-2*delta,ODn_3+2*b+2*w),180e-3,deg2rad(2),le_min];         % lower bound
%ubd = [min(a,D-a/sin(pi/n)),min(a,D-a/sin(pi/n))-2*w-2*delta,min(a,D-a/sin(pi/n))-4*w-4*delta,200e-3,deg2rad(15),10*le_min]; % upper bound
%ubd = [4.95e-3,4.1e-3,200e-3,deg2rad(7.5),203e-3]; % upper bound
%lbd = [2.38e-3,1.98e-3,192e-3,deg2rad(3.7),le_min];         % lower bound
lbd = [2.8e-3,2.4e-3,1.9e-3,195e-3,deg2rad(2),40e-3];         % lower bound
ubd = [5.0e-3,4.3e-3,3.2e-3,205e-3,deg2rad(9),215e-3]; % upper bound
nonlcon = @nonlinear_constraints_2;               %nonlinear constraints
%PopulationSize_Data= 200000;                      %gamultiobj setting parameters
PopulationSize_Data= 40000;                      %gamultiobj setting parameters
FunctionTolerance_Data = 1e-9;                  %gamultiobj setting parameters
ConstraintTolerance_Data = 1e-9;                %gamultiobj setting parameters

options = optimoptions('gamultiobj');           %default options
%modify options setting
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'FunctionTolerance', FunctionTolerance_Data);
options = optimoptions(options,'ConstraintTolerance', ConstraintTolerance_Data);
options = optimoptions(options,'CreationFcn', @gacreationnonlinearfeasible);
options = optimoptions(options,'CrossoverFcn', {  @crossoverintermediate [] });
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'PlotFcn', { @gaplotpareto });   % plot pareto front

% Using parallel computing
options = optimoptions(options,'UseParallel',true);

% run galmultiobj function and get the pareto front
[x_gal,fval_gal] = gamultiobj(func,nvars,Aineq,bineq,Aeq,beq,lbd,ubd,nonlcon,options);
toc

% !! You can click the "stop" button while Pareto Front appears
