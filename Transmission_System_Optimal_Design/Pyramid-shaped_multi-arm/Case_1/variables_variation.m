%% Variation of Optimization Varables
% By Chao Zhang
% Date: 2022/11/20

clear;clc;


%% Load CTR params
CTR_params;


%% Optimization variables
ODs_1 = optimvar('ODs_1','LowerBound',ODn_1+2*b+2*w,'UpperBound',min(a,D-a/sin(pi/n)));
ODs_2 = optimvar('ODs_2','LowerBound',max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),'UpperBound',min(a,D-a/sin(pi/n))-2*w-2*delta);
lb = optimvar('lb','LowerBound',50e-3,'UpperBound',200e-3);
alpha = optimvar('alpha','LowerBound',deg2rad(3),'UpperBound',deg2rad(15));
le = optimvar('le','LowerBound',le_min,'UpperBound',10*le_min);
beta_1 = optimvar('beta_1','LowerBound',0,'UpperBound',pi);
beta_2 = optimvar('beta_2','LowerBound',0,'UpperBound',pi); 


%% Geometric Constraints - Relations of tube diameters
IDs_2 = ODn_2+2*b;                      % % Inner diameter of SS tube 2 (inner tube)
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


%% Calculate optimization variables looply
% Giving a varying zeta_max, solving min(le+lb)
X = 100;                                    % Calculate X times
zeta_init = 0.058677334853678;              % Initial zeta_max, rad (change it manually)
zeta_step = (zeta_max - zeta_init)/(X-1);
zeta_max = zeros(X,1);
sol_data = zeros(X,7);

for i=1:X
    zeta_max(i) = zeta_init + (i-1)*zeta_step;
    
    prob = optimproblem;

    % Optimization function
    prob.Objective = le+lb;    % Compactness
    prob.ObjectiveSense = 'minimize';

    % Design constraints
    cons_1 = ODs_2<=ODs_1-2*w-2*delta;  % Tube diameters
    cons_2 = epsilon<=epsilon_e;        % Elastic strain limit
    cons_3 = zeta_1==zeta_2;            % Maximum relative twist angle
    cons_4 = zeta_2<=zeta_max(i);       % Maximum relative twist angle
    cons_5 = d1>=0;                     % collision avoidance
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

    % Record results
    sol_data(i,1) = sol.ODs_1;
    sol_data(i,2) = sol.ODs_2;
    sol_data(i,3) = sol.lb;
    sol_data(i,4) = sol.alpha;
    sol_data(i,5) = sol.le;
    sol_data(i,6) = sol.beta_1;
    sol_data(i,7) = sol.beta_2;
end


%% Plot the variation of the optimization variables
% zeta vs ODs_1 and ODs_2 
subplot(2,2,1)
maker_idx = [1,10,20,30,40,50,60,70,80,90,100];
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,1),'-s','MarkerIndices',maker_idx,'LineWidth',1.5)
hold on
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,2),'-^','MarkerIndices',maker_idx,'LineWidth',1.5)
xlim([3 8]);
ylim([3.4 5]);
set(gca,'YTick',[3.4:0.4:5]);
grid on;
legend('ODs_1','ODs_2');
xlabel('F_1(x) (°)');
ylabel('Diameter (mm)');


% zeta vs lb and le 
subplot(2,2,2)
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,3),'-s','MarkerIndices',maker_idx,'LineWidth',1.5);
hold on
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,5),'-^','MarkerIndices',maker_idx,'LineWidth',1.5);
xlim([3 8]);
ylim([80 240]);
set(gca,'YTick',[80:40:240]);
grid on;
legend('l_b','l_e');
xlabel('F_1(x) (°)');
ylabel('Length (mm)');


% zeta vs alpha
subplot(2,2,3)
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,4)),'-s','MarkerIndices',maker_idx,'LineWidth',1.5)
xlim([3 8]);
ylim([3.6 4.6]);
set(gca,'YTick',[3.6:0.2:4.6]);
grid on;
xlabel('F_1(x) (°)');
ylabel('\alpha (°)');


% zeta vs beta
subplot(2,2,4)
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,6)),'-s','MarkerIndices',maker_idx,'LineWidth',1.5)
hold on
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,7)),'-^','MarkerIndices',maker_idx+2,'LineWidth',1.5)
xlim([3 8]);
ylim([0 160]);
set(gca,'YTick',[0:40:160]);
grid on;
legend('\beta_1','\beta_2');
xlabel('F_1(x) (°)');
ylabel('\beta_i (°)');

%% 导出设置
%大小 - 20x14cm
%渲染 - 600dpi
%字体 - 14，Times New Roman
