%% Optimal Design Framework for Pyramid-shaped Transmission System
% This case demonstrates the design optimization of a triple-arm CTR with 2 tubes (balanced pair) in each arm.
% 1. Establish the design constraints
% 2. Formulate optimization problem
% 3. Getting the Pareto front
% by Chao Zhang
% Date：2022/12/8

clear;clc;
syms ODs_1 ODs_2 lb alpha le beta_1 beta_2 real    % % Design variables


%% Load CTR params
CTR_params;


%% Geometric Constraints - Relations of tube diameters
IDs_2 = ODn_2+2*b;                      % Inner diameter of SS tube 2 (inner tube)
%ODs_1 <= min(a,D-a/cos(pi/2-pi/n));
%ODs_1 >= ODn_1+2*b+2*w;
%ODs_2 >= max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w)
%ODs_2 <= ODs_1-2*w-2*delta


%% Geometric Constraints - Shape and dimensions of SS tubes
v = 2/3 * lb * alpha;                   % Maximum deflection of SS tube
L1 = delta_L_1+le+lb+ls;                % Length of SS tube 1
L2 = delta_L_2+le+lb+ls-la;             % Length of SS tube 2
% le >= le_min


%% Geometric Constraints - Collision avoidance of actuation units
% Coordinates of tube points
O = [0;0;0];
O1  = [sqrt(3)/3*a;0;0];
O2  = [-sqrt(3)/6*a;a/2;0];
O3  = [-sqrt(3)/6*a;-a/2;0];

A = [0;0;ls];
A1 = [sqrt(3)/3*a;0;ls];
A2 = [-sqrt(3)/6*a;a/2;ls];
A3 = [-sqrt(3)/6*a;-a/2;ls];

B = [0;0;ls+lb];
B1 = [sqrt(3)/3*a+v;0;ls+lb];
B2 = [-sqrt(3)/6*a+v*cos(beta_1);a/2+v*sin(beta_1);ls+lb];
B3 = [-sqrt(3)/6*a+v*cos(beta_1+beta_2);-a/2+v*sin(beta_1+beta_2);ls+lb];

C = [0;0;ls+lb+le*cos(alpha)];
C1 = [sqrt(3)/3*a+v+le*sin(alpha);0;ls+lb+le*cos(alpha)];
C2 = [-sqrt(3)/6*a+(v+le*sin(alpha))*cos(beta_1);a/2+(v+le*sin(alpha))*sin(beta_1);ls+lb+le*cos(alpha)];
C3 = [-sqrt(3)/6*a+(v+le*sin(alpha))*cos(beta_1+beta_2);-a/2+(v+le*sin(alpha))*sin(beta_1+beta_2);ls+lb+le*cos(alpha)];

% line vector
C1C2 = C2-C1;
C2C3 = C3-C2;
C3C1 = C1-C3;
d_C1C2 = sqrt(C1C2(1)^2+C1C2(2)^2);
d_C2C3 = sqrt(C2C3(1)^2+C2C3(2)^2);
d_C3C1 = sqrt(C3C1(1)^2+C3C1(2)^2);

% Normal vector of the plane \kappa_j
n1 = (C1-B1)/le;
n2 = (C2-B2)/le;
n3 = (C3-B3)/le;

% Included angle between the line and the plane
sin_gama12 = dot(n1,C1C2)/d_C1C2;
sin_gama21 = dot(n2,-C1C2)/d_C1C2;
sin_gama23 = dot(n2,C2C3)/d_C2C3;
sin_gama32 = dot(n3,-C2C3)/d_C2C3;
sin_gama31 = dot(n3,C3C1)/d_C3C1;
sin_gama13 = dot(n1,-C3C1)/d_C3C1;

% Projecting the radius ra
d_C1D12 = ra*sqrt(1-sin_gama12^2);
d_C2D21 = ra*sqrt(1-sin_gama21^2);
d_C2D23 = ra*sqrt(1-sin_gama23^2);
d_C3D32 = ra*sqrt(1-sin_gama32^2);
d_C3D31 = ra*sqrt(1-sin_gama31^2);
d_C1D13 = ra*sqrt(1-sin_gama13^2);

% Minimum distance between adjacent actuation units
d1 = ( d_C1C2 - d_C1D12 - d_C2D21 );
d2 = ( d_C2C3 - d_C2D23 - d_C3D32 );
d3 = ( d_C3C1 - d_C3D31 - d_C1D13 );


%% Deformation Constraints - Elastic strain limit of bending
epsilon = ODs_1*alpha/lb;               % Maximum bendign strain in SS tube
%epsilon <= epsilon_e


%% Deformation Constraints - Torsional twisting of SS tubes
% Inner SS tube
J2 = pi*(ODs_2^4-IDs_2^4)/32;           % Cross-sectional polar moment of inertia, m^4
kz_2 = Gs*J2;                           % Torsion stiffness, m^4*Pa
zeta_2 = tau_2_max*(L2)/kz_2;           % Maximum relative twist angle, rad

% Outer SS tube
IDs_1 = ODs_2+2*delta;                  % Inner diameter of outter SS tube
J1 = pi*(ODs_1^4-IDs_1^4)/32;           % Cross-sectional polar moment of inertia, m^4
kz_1 = Gs*J1;                           % Torsion stiffness, m^4*Pa
zeta_1 = tau_1_max*(L1)/kz_1;           % Maximum relative twist angle, rad


%% function and constraints generation
% !! Update "multiobjective_func.m" and "nonlinear_constraints.m" manually

% F(1) = zeta_i
Func_1 = (180/pi)*zeta_1

% F(2) = le+lb
Func_2 = le + lb + ls - la + delta_L_2

% ODs_2 <= ODs_1-2*w-2*delta
Cons_1 = ODs_2 - ODs_1 + 2*w + 2*delta

% epsilon <= epsilon_e
Cons_2 = epsilon - epsilon_e

% zeta_i<=zeta_max
Cons_3 = (180/pi)*(zeta_2 - zeta_max)

% d1>=0 -> -d1<=0 (make beta_1 = 2*pi/3)
Cons_4 = - subs(d1,beta_1,2*pi/n)

% zeta_1=zeta_2
Cons_5 = (180/pi)*(zeta_1 - zeta_2)


%% Formulation of the multiobjective optimization problem
tic
func = @multiobjective_func;            %objective function
nvars = 5;                              %numbers of optimization variables
Aineq = [-1,1,0,0,0];                   %linear constraint: Ax<=b
bineq = -2*w-2*delta;                   %linear inequally constraint: Ax<=b
Aeq = [];                               %linear equally constraints:Aeq x = beq
beq = [];
lbd = [ODn_1+2*b+2*w,max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),50e-3,deg2rad(2),le_min];         % lower bound
ubd = [min(a,D-a/sin(pi/n)),min(a,D-a/sin(pi/n))-2*w-2*delta,200e-3,deg2rad(10),10*le_min]; % upper bound
nonlcon = @nonlinear_constraints;               %nonlinear constraints
PopulationSize_Data= 80000;                      %gamultiobj setting parameters
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
