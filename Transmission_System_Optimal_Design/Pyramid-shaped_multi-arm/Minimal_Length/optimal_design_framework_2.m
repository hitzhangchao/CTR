%% Optimal Design Framework for Pyramid-shaped Transmission System
% 主优化程序
% by Zhang Chao
% Date：2022/10/20

clear;clc;
syms ODs_1 ODs_2 lb alpha le beta_1 beta_2 real    % Design variables: ODs_i-SS管外径；le-近端直线长度；lb-圆弧长度；alpha-弯曲角度

%% Load CTR params
CTR_params_2;

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


%% function and constraints generation - 多目标优化部分需半自动手敲，以提升运行速度
% F(1) = zeta_i
Func_1 = (180/pi)*zeta_1

% F(2) = le+lb
Func_2 = le + lb

% ODs_2 <= ODs_1-2*w-2*delta
Cons_1 = ODs_2 - ODs_1 + 2*w + 2*delta

% epsilon <= epsilon_e
Cons_2 = epsilon - epsilon_e

% zeta_i<=zeta_max
Cons_3 = (180/pi)*(zeta_2 - zeta_max)

% d1>=0 -> -d1<=0 (make beta_1 = pi/2)
Cons_4 = - subs(d1,beta_1,2*pi/n)

% zeta_1=zeta_2
Cons_5 = (180/pi)*(zeta_1 - zeta_2)


%% Formulation of the multiobjective optimization problem
tic
func = @multiobjective_func_2;            %objective function
nvars = 5;                              %numbers of optimization variables
Aineq = [-1,1,0,0,0];                   %linear constraint: Ax<=b
bineq = -2*w-2*delta;                   %linear inequally constraint: Ax<=b
Aeq = [];                               %linear equally constraints:Aeq x = beq
beq = [];
%lbd = [ODn_1+2*b+2*w,max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),178e-3,deg2rad(5),le_min];         % lower bound
%ubd = [min(a,D-a/sin(pi/n)),min(a,D-a/sin(pi/n))-2*w-2*delta,200e-3,deg2rad(10),4*le_min]; % upper bound
%ubd = [4.95e-3,4.1e-3,200e-3,deg2rad(7.5),203e-3]; % upper bound
%lbd = [2.38e-3,1.98e-3,192e-3,deg2rad(3.7),le_min];         % lower bound
lbd = [2.3e-3,1.90e-3,50e-3,deg2rad(3),le_min];         % lower bound
ubd = [5.0e-3,4.2e-3,200e-3,deg2rad(8),220e-3]; % upper bound
nonlcon = @nonlinear_constraints_2;               %nonlinear constraints
PopulationSize_Data= 200000;                      %gamultiobj setting parameters
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

% 并行加速
options = optimoptions(options,'UseParallel',true);

% run galmultiobj function and get the pareto front
[x_gal,fval_gal] = gamultiobj(func,nvars,Aineq,bineq,Aeq,beq,lbd,ubd,nonlcon,options);
toc

%% Pareto front图像输出
% 这曲线质量比之前高多了，运行速度也快多了！
% 输出高12mm，宽20mm