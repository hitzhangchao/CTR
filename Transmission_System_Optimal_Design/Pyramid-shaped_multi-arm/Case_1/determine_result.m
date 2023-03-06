%% Determine the Final Optimizing Result
% by Chao Zhang
% Date：2022/11/21

clear;clc;

%% Load CTR params
CTR_params;

%% Choose a Preferred Point on Pareto front in which
zeta_max = deg2rad(4);


%% Optimization Variables
ODs_1 = optimvar('ODs_1','LowerBound',ODn_1+2*b+2*w,'UpperBound',min(a,D-a/sin(pi/n)));
ODs_2 = optimvar('ODs_2','LowerBound',max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),'UpperBound',min(a,D-a/sin(pi/n))-2*w-2*delta);
lb = optimvar('lb','LowerBound',50e-3,'UpperBound',200e-3);
alpha = optimvar('alpha','LowerBound',deg2rad(3),'UpperBound',deg2rad(15));
le = optimvar('le','LowerBound',le_min,'UpperBound',10*le_min);
beta_1 = optimvar('beta_1','LowerBound',0,'UpperBound',pi);
beta_2 = optimvar('beta_2','LowerBound',0,'UpperBound',pi); 


%% Geometric Constraints - Relations of tube diameters
IDs_2 = ODn_2+2*b;                      % Inner diameter of SS tube 2 (inner tube)
%ODs_1 <= min(a,D-a/cos(pi/2-pi/n));
%ODs_1 >= ODn_1+2*b+2*w;
%ODs_2 >= max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w)
%ODs_2 <= ODs_1-2*w-2*delta


%% Geometric Constraints - Shape and length of SS tubes
v = 2/3 * lb * alpha;                   % Maximum deflection of SS tube
L1 = delta_L_1+le+lb+lc;                % Length of SS tube 1
L2 = delta_L_2+le+lb+lc;                % Length of SS tube 2
% le >= le_min


%% Geometric Constraints - Collision avoidance of actuation units
% Coordinates of tube points
O = [0;0;0];
O1  = [sqrt(3)/3*a;0;0];
O2  = [-sqrt(3)/6*a;a/2;0];
O3  = [-sqrt(3)/6*a;-a/2;0];

A = [0;0;lc];
A1 = [sqrt(3)/3*a;0;lc];
A2 = [-sqrt(3)/6*a;a/2;lc];
A3 = [-sqrt(3)/6*a;-a/2;lc];

B = [0;0;lc+lb];
B1 = [sqrt(3)/3*a+v;0;lc+lb];
B2 = [-sqrt(3)/6*a+v*cos(beta_1);a/2+v*sin(beta_1);lc+lb];
B3 = [-sqrt(3)/6*a+v*cos(beta_1+beta_2);-a/2+v*sin(beta_1+beta_2);lc+lb];

C = [0;0;lc+lb+le*cos(alpha)];
C1 = [sqrt(3)/3*a+v+le*sin(alpha);0;lc+lb+le*cos(alpha)];
C2 = [-sqrt(3)/6*a+(v+le*sin(alpha))*cos(beta_1);a/2+(v+le*sin(alpha))*sin(beta_1);lc+lb+le*cos(alpha)];
C3 = [-sqrt(3)/6*a+(v+le*sin(alpha))*cos(beta_1+beta_2);-a/2+(v+le*sin(alpha))*sin(beta_1+beta_2);lc+lb+le*cos(alpha)];

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
%di>=0


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
J1 = pi*(ODs_1^4-IDs_1^4)/32;           % Cross-sectional polar moment of inertia
kz_1 = Gs*J1;                           % Torsion stiffness
zeta_1 = tau_1_max*(L1)/kz_1;           % Maximum relative twist angle


%% Formulating the Optimization Problem
prob = optimproblem;

% Optimization function
choose_func = 2;
if choose_func == 1
    prob.Objective = zeta_1;    % Maxium relative twist angle
else
    prob.Objective = le+lb;     % System compactness
end
prob.ObjectiveSense = 'minimize';

% Design constraints
cons_1 = ODs_2<=ODs_1-2*w-2*delta;  % Tube diameters
cons_2 = epsilon<=epsilon_e;        % Elastic strain limit
cons_3 = zeta_1==zeta_2;            % Maximum relative twist angle
cons_4 = zeta_2<=zeta_max;          % Maximum relative twist angle
cons_5 = d1>=0;                     % Collision avoidance
cons_6 = d2>=0;
cons_7 = d3>=0;

prob.Constraints.cons1 = cons_1;
prob.Constraints.cons2 = cons_2;
prob.Constraints.cons3 = cons_3;
prob.Constraints.cons4 = cons_4;
prob.Constraints.cons5 = cons_5;
prob.Constraints.cons6 = cons_6;
prob.Constraints.cons7 = cons_7;

show(prob)

% Initial values
initialpt.ODs_1 = ODn_1;
initialpt.ODs_2 = ODn_2;
initialpt.lb = 0;
initialpt.alpha = 0;
initialpt.le = le_min;
initialpt.beta_1 = 2*pi/n;
initialpt.beta_2 = 2*pi/n;

% Solve the problem
[sol,fval] = solve(prob,initialpt)

% Check the result
nx = norm([sol.ODs_1 sol.ODs_2 sol.lb sol.le sol.alpha sol.beta_1 sol.beta_2])

%% Optimal Design Results
ODs_1 = sol.ODs_1
ODs_2 = sol.ODs_2
lb = sol.lb
alpha_deg =rad2deg(sol.alpha)
le = sol.le
beta_1 = rad2deg(sol.beta_1)
beta_2 = rad2deg(sol.beta_2)

IDs_1 = ODs_2+2*delta
IDs_2
L1 = delta_L_1+le+lb+lc
L2 = delta_L_2+le+lb+lc
zeta_1_deg = rad2deg(evaluate(zeta_1,sol))
zeta_2_deg = rad2deg(evaluate(zeta_2,sol))

