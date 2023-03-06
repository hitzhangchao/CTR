%%  The varing of the design varables
%   Coding by: Zhang Chao
%   Date: 2022/11/20
%   This case demonstrates the design optimization of a 4 arm CTR with 2 tubes in each arm.

clear;clc;

%% Load CTR params
CTR_params_2;


%% Optimization variables
ODs_1 = optimvar('ODs_1','LowerBound',ODn_1+2*b+2*w,'UpperBound',min(a,D-a/sin(pi/n)));
ODs_2 = optimvar('ODs_2','LowerBound',max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w),'UpperBound',min(a,D-a/sin(pi/n))-2*w-2*delta);
lb = optimvar('lb','LowerBound',50e-3,'UpperBound',200e-3);
alpha = optimvar('alpha','LowerBound',deg2rad(3),'UpperBound',deg2rad(10));
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


%% Calculate optimization variables looply
%给定不同的zeta_max值，以le+lb为优化目标函数，循环计算 X 次
X = 100;
zeta_init = 0.002914194585748;             % 扭转角初始最小值
zeta_step = (zeta_max - zeta_init)/(X-1);   % zeta_max每次循环增加的步长
zeta_max = zeros(X,1);                      %为数组预分配内存
sol_data = zeros(X,8);                      %为记录优化数据的数组分配内存

for i=1:X
    zeta_max(i) = zeta_init + (i-1)*zeta_step;
    
    prob = optimproblem;
    %目标函数
    prob.Objective = le+lb;    % 尺寸最compact
    prob.ObjectiveSense = 'minimize';

    %约束条件
    cons_1 = ODs_2<=ODs_1-2*w-2*delta;  % SS管外径ODs_2的尺寸约束
    cons_2 = epsilon<=epsilon_e;        % 线性应变约束
    cons_3 = zeta_1==zeta_2;            % 各管相对扭转角相等
    cons_4 = zeta_2<=zeta_max(i);          % 相对扭转角小于最大允许值
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

    %记录每次优化的结果
    sol_data(i,1) = sol.ODs_1;
    sol_data(i,2) = sol.ODs_2;
    sol_data(i,3) = sol.lb;
    sol_data(i,4) = sol.alpha;
    sol_data(i,5) = sol.le;
    sol_data(i,6) = pi/2;
    sol_data(i,7) = pi/2;
    sol_data(i,8) = pi/2;
end


%% 绘制关系图线
%% zeta - ODs_1, ODs_1 
subplot(2,2,1)
maker_idx = [1,10,20,30,40,50,60,70,80,90,100];
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,1),'-s','MarkerIndices',maker_idx,'LineWidth',1.5)
hold on
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,2),'-^','MarkerIndices',maker_idx,'LineWidth',1.5)
xlim([0 3]);
ylim([1.8 5]);
set(gca,'YTick',[1.8:0.8:5]);
grid on;
legend('ODs_1','ODs_2');
xlabel('F_1(x) (°)');
ylabel('Diameter (mm)');

%% zeta - lb, le 
subplot(2,2,2)
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,3),'-s','MarkerIndices',maker_idx,'LineWidth',1.5);
hold on
plot(rad2deg(zeta_max(:,1)),1000*sol_data(:,5),'-^','MarkerIndices',maker_idx,'LineWidth',1.5);
xlim([0 3]);
ylim([0 240]);
set(gca,'YTick',[0:60:240]);
grid on;
legend('l_b','l_e');
xlabel('F_1(x) (°)');
ylabel('Length (mm)');

%% zeta - alpha
subplot(2,2,3)
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,4)),'-s','MarkerIndices',maker_idx,'LineWidth',1.5)
xlim([0 3]);
ylim([3.5 7.5]);
set(gca,'YTick',[3.5:1:7.5]);
grid on;
xlabel('F_1(x) (°)');
ylabel('\alpha (°)');

%% zeta - beta 
subplot(2,2,4)
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,6)),'-s','MarkerIndices',maker_idx,'LineWidth',1.5)
hold on
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,7)),'-^','MarkerIndices',maker_idx+2,'LineWidth',1.5)
hold on
plot(rad2deg(zeta_max(:,1)),rad2deg(sol_data(:,8)),'-o','MarkerIndices',maker_idx+4,'LineWidth',1.5)
xlim([0 3]);
ylim([0 120]);
set(gca,'YTick',[0:30:120]);
grid on;
legend('\beta_1','\beta_2','\beta_3');
xlabel('F_1(x) (°)');
ylabel('\beta_i (°)');

%% 导出设置
%大小 - 20x14cm
%渲染 - 600dpi
%字体 - 14，Times New Roman


