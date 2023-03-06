%% Determine the Final Optimizing Result
% 根据Pareto front，选择最优点，确定各设计参数
% by Zhang Chao
% Date：2022/11/20

clear;clc;

%% Load CTR params
CTR_params_2;


%% Optimization variables
ODs_1 = optimvar('ODs_1','LowerBound',ODn_1+2*b+2*w,'UpperBound',min(a,D-a/sin(pi/n)));
ODs_2 = optimvar('ODs_2','LowerBound',max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),'UpperBound',min(a,D-a/sin(pi/n))-2*w-2*delta);
lb = optimvar('lb','LowerBound',50e-3,'UpperBound',200e-3);
alpha = optimvar('alpha','LowerBound',deg2rad(3),'UpperBound',deg2rad(15));
le = optimvar('le','LowerBound',le_min,'UpperBound',10*le_min);


%% Geometric Constraints - Relations of tube diameters
IDs_2 = ODn_2+2*b;                      % 最内SS管内径
%ODs_1 <= min(a,D-a/cos(pi/2-pi/n));
%ODs_1 >= ODn_1+2*b+2*w;
%ODs_2 >= max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w)
%ODs_2 <= ODs_1-2*w-2*delta


%% Geometric Constraints - Shape and dimensions of SS tubes
v = 2/3 * lb * alpha;                   % 最大挠曲变形
L1 = delta_L_1+le+lb+lc;                % 管1总长度
L2 = delta_L_2+le+lb+lc;                % 管2总长度
% le >= le_min


%% Geometric Constraints - Collision avoidance of actuation units
% 建立坐标系，空间各点坐标 - 这里直接赋值beta=90°
O = [0;0;0];
O1 = [sqrt(2)/2*a;0;0];
O2 = [0;sqrt(2)/2*a;0];
O3 = [-sqrt(2)/2*a;0;0];
O4 = [0;-sqrt(2)/2*a;0];       % 4 arm，需要增加一个

A = [0;0;lc];
A1 = [sqrt(2)/2*a;0;lc];
A2 = [0;sqrt(2)/2*a;lc];
A3 = [-sqrt(2)/2*a;0;lc];
A4 = [0;-sqrt(2)/2*a;lc];

B = [0;0;lc+lb];
B1 = [sqrt(2)/2*a+v;0;lc+lb];
B2 = [0;sqrt(2)/2*a+v;lc+lb];
B3 = [-sqrt(2)/2*a-v;0;lc+lb];
B4 = [0;-sqrt(2)/2*a-v;lc+lb];

C = [0;0;lc+lb+le*cos(alpha)];
C1 = [sqrt(2)/2*a+v+le*sin(alpha);0;lc+lb+le*cos(alpha)];
C2 = [0;sqrt(2)/2*a+v+le*sin(alpha);lc+lb+le*cos(alpha)];
C3 = [-sqrt(2)/2*a-v-le*sin(alpha);0;lc+lb+le*cos(alpha)];
C4 = [0;-sqrt(2)/2*a-v-le*sin(alpha);lc+lb+le*cos(alpha)];

% actuation module中心点连线矢量
C1C2 = C2-C1;
C2C3 = C3-C2;
C3C4 = C4-C3;
C4C1 = C1-C4;
d_C1C2 = sqrt(C1C2(1)^2+C1C2(2)^2);
d_C2C3 = sqrt(C2C3(1)^2+C2C3(2)^2);
d_C3C4 = sqrt(C3C4(1)^2+C3C4(2)^2);
d_C4C1 = sqrt(C4C1(1)^2+C4C1(2)^2);

% actuation module平面法向量
n1 = (C1-B1)/le;
n2 = (C2-B2)/le;
n3 = (C3-B3)/le;
n4 = (C4-B4)/le;

% 连线矢量和平面夹角
sin_gama12 = dot(n1,C1C2)/d_C1C2;
sin_gama21 = dot(n2,-C1C2)/d_C1C2;
sin_gama23 = dot(n2,C2C3)/d_C2C3;
sin_gama32 = dot(n3,-C2C3)/d_C2C3;
sin_gama34 = dot(n3,C3C4)/d_C3C4;
sin_gama43 = dot(n4,-C3C4)/d_C3C4;
sin_gama41 = dot(n4,C4C1)/d_C4C1;
sin_gama14 = dot(n1,-C4C1)/d_C4C1;

% 圆面线段投影
d_C1D12 = ra*sqrt(1-sin_gama12^2);
d_C2D21 = ra*sqrt(1-sin_gama21^2);
d_C2D23 = ra*sqrt(1-sin_gama23^2);
d_C3D32 = ra*sqrt(1-sin_gama32^2);
d_C3D34 = ra*sqrt(1-sin_gama34^2);
d_C4D43 = ra*sqrt(1-sin_gama43^2);
d_C4D41 = ra*sqrt(1-sin_gama41^2);
d_C1D14 = ra*sqrt(1-sin_gama14^2);

d1 = d_C1C2 - d_C1D12 - d_C2D21;
d2 = d_C2C3 - d_C2D23 - d_C3D32;
d3 = d_C3C4 - d_C3D34 - d_C4D43;
d4 = d_C4C1 - d_C4D41 - d_C1D14;


%% Deformation Constraints - Elastic strain limit of bending
epsilon = ODs_1*alpha/lb;               % SS管应变
%epsilon <= epsilon_e                    %线性应变约束


%% Deformation Constraints - Torsional twisting of SS tubes
% 内管
J2 = pi*(ODs_2^4-IDs_2^4)/32;           % 极惯性矩，m^4
kz_2 = Gs*J2;                           % 扭转刚度，m^4*Pa
zeta_2 = tau_2_max*(L2)/kz_2;           % 扭转角度,rad

% 外管
IDs_1 = ODs_2+2*delta;                  % 外管内径
J1 = pi*(ODs_1^4-IDs_1^4)/32;           % 极惯性矩，m^4
kz_1 = Gs*J1;                           % 扭转刚度，m^4*Pa
zeta_1 = tau_1_max*(L1)/kz_1;           % 扭转角度,rad


%% 创建优化问题
prob = optimproblem;

%目标函数
choose_func = 1;
if choose_func == 1
    prob.Objective = zeta_1;    % 相对扭转角最小
else
    prob.Objective = le+lb;    % 尺寸最compact
end
prob.ObjectiveSense = 'minimize';

%约束条件
cons_1 = ODs_2<=ODs_1-2*w-2*delta;  % SS管外径ODs_2的尺寸约束
cons_2 = epsilon<=epsilon_e;        % 线性应变约束
cons_3 = zeta_1==zeta_2;            % 各管相对扭转角相等
cons_4 = zeta_2<=zeta_max;          % 相对扭转角小于最大允许值
cons_5 = d1>=0;                      %collision avoidance
cons_6 = d2>=0;
cons_7 = d3>=0;
cons_8 = d4>=0;

%在问题中包含约束 - 其余的约束直接在定义优化变量的时候在上下界处给出了
prob.Constraints.cons1 = cons_1;
prob.Constraints.cons2 = cons_2;
prob.Constraints.cons3 = cons_3;
prob.Constraints.cons4 = cons_4;
prob.Constraints.cons5 = cons_5;
prob.Constraints.cons6 = cons_6;
prob.Constraints.cons7 = cons_7;
prob.Constraints.cons8 = cons_8;

%检查此优化问题
show(prob)

%初值
initialpt.ODs_1 = ODn_1;
initialpt.ODs_2 = ODn_2;
initialpt.lb = 0;
initialpt.alpha = 0;
initialpt.le = le_min;

%求解问题
[sol,fval] = solve(prob,initialpt)

%检查解
nx = norm([sol.ODs_1 sol.ODs_2 sol.lb sol.le sol.alpha])                        %范数要确保小于或等于1

%% Optimal Design Results
ODs_1 = sol.ODs_1
ODs_2 = sol.ODs_2
lb = sol.lb
alpha_deg =rad2deg(sol.alpha)
le = sol.le

IDs_1 = ODs_2+2*delta
IDs_2
L1 = delta_L_1+le+lb+lc
L2 = delta_L_2+le+lb+lc
zeta_1_deg = rad2deg(evaluate(zeta_1,sol))
zeta_2_deg = rad2deg(evaluate(zeta_2,sol))

