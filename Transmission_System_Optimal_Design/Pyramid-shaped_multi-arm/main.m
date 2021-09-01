%% 多臂CTR平台传动优化设计 - 以3臂，每臂为一对 balance pair为例
% 优化问题中，我们让每个SS管的最大扭转角相等，并以扭转角最小化为目标函数，这样就能不使用多目标优化

clc;clear;
syms ODs_1 ODs_2 le lb alpha beta_1 beta_2 real  %ODs_i-SS管外径；le-近端直线长度；lb-圆弧长度；alpha-弯曲角度

%% NiTi设计需求 - 双臂3管为例 
ODn_2 = 0.0025;                         %NiTi管外径,m
ODn_1 = 0.0029;                         %最外NiTi管外径,m
a = 0.0072;                             %管的中心距（以双臂为例）,m
b = 0.0001;                             %SS管和NiTi管最小粘接间隙,m
w = 0.0001;                             %SS管最小壁厚
delta = 0.0001;                         %SS管间最小游动间隙
delta_L_1 = 0.050;                  	% SS管1电机端到lp起点的距离
delta_L_2 = 0.110;                      % SS管2电机端到le起点的距离
tau_2 = 0.2;                         	%需传递力矩，Nm
tau_1 = 0.22;                         	%需传递力矩，Nm
lc = 0.080;                             %SS即管伸出直线段的长度,准直器
IDs_2 = ODn_2+2*b;                      %最内SS管内径

% 管尺寸约束
%ODn_1+2*b+2*w <= ODs_1 <= a
%max(ODn_1+2*b-2*delta,ODn_2+2*b+2*w) <= ODs_2 <=ODs_1-2*w-2*delta
%max(ODn_2+2*b-2*delta,ODn_3+2*b+2*w) <= ODs_3 <= ODs_2-2*w-2*delta


%% 假设ss管弯曲段近似为圆弧
r = ODs_1/2;                            %最外SS管半径
R = lb/alpha;                           %圆弧段弯曲半径
epsilon = r/R;                          %SS管应变


%% 管线性应变约束
epsilon_max = 0.0011;                    %SS管线性应变段最大允许应变（如果超过，管子就可能塑性变形）
%约束条件2
% epsilon <= epsilon_max                %线性应变约束


%% 输入直线段L_s的长度约束
le_min = 0.050;                         %lp最小长度
%约束条件3
% le >= le_min                          


%% 驱动端不干涉约束
lm = 0.012;                             %驱动端最大宽度，m
d_min = 0.010;                          %驱动端避免碰撞的最小距离
% 建立坐标系，图中各点坐标
%A = [0;a/2];
%B = [-R*sin(alpha);a/2+R*(1-cos(alpha))];
%C = [-R*sin(alpha)-le*cos(alpha);a/2+R*(1-cos(alpha))+le*sin(alpha)];
%D = [-R*sin(alpha)-le*cos(alpha)-lm*sin(alpha);a/2+R*(1-cos(alpha))+le*sin(alpha)-lm*cos(alpha)];
%Dy = 2*D(2);                            %驱动端两最近点距离
% Dy >= d_min                               %约束条件4

% test - 距离计算是否正确？
%beta_1 = 2/3*pi;
%beta_2 = 2/3*pi;
%lb = 0.200;
%alpha = deg2rad(10);
%le = 0.050;

%空间各点坐标
O = [0;0;0];
O1  = [sqrt(3)/3*a;0;0];
O2  = [-sqrt(3)/6*a;a/2;0];
O3  = [-sqrt(3)/6*a;-a/2;0];

A = [0;0;lc];
A1 = [sqrt(3)/3*a;0;lc];
A2 = [-sqrt(3)/6*a;a/2;lc];
A3 = [-sqrt(3)/6*a;-a/2;lc];

B = [0;0;lc+R*sin(alpha)];
B1 = [sqrt(3)/3*a+R*(1-cos(alpha));0;lc+R*sin(alpha)];
B2 = [-sqrt(3)/6*a+R*(1-cos(alpha))*cos(beta_1);a/2+R*(1-cos(alpha))*sin(beta_1);lc+R*sin(alpha)];
B3 = [-sqrt(3)/6*a+R*(1-cos(alpha))*cos(beta_1+beta_2);-a/2+R*(1-cos(alpha))*sin(beta_1+beta_2);lc+R*sin(alpha)];

C = [0;0;lc+R*sin(alpha)+le*cos(alpha)];
C1 = [sqrt(3)/3*a+R*(1-cos(alpha))+le*sin(alpha);0;lc+R*sin(alpha)+le*cos(alpha)];
C2 = [-sqrt(3)/6*a+(R*(1-cos(alpha))+le*sin(alpha))*cos(beta_1);a/2+(R*(1-cos(alpha))+le*sin(alpha))*sin(beta_1);lc+R*sin(alpha)+le*cos(alpha)];
C3 = [-sqrt(3)/6*a+(R*(1-cos(alpha))+le*sin(alpha))*cos(beta_1+beta_2);-a/2+(R*(1-cos(alpha))+le*sin(alpha))*sin(beta_1+beta_2);lc+R*sin(alpha)+le*cos(alpha)];

% actuation module中心点连线矢量
C1C2 = C2-C1;
C2C3 = C3-C2;
C3C1 = C1-C3;
d_C1C2 = sqrt(C1C2(1)^2+C1C2(2)^2);
d_C2C3 = sqrt(C2C3(1)^2+C2C3(2)^2);
d_C3C1 = sqrt(C3C1(1)^2+C3C1(2)^2);

% actuation module平面法向量
n1 = simplify(C1-B1)/le;
n2 = simplify(C2-B2)/le;
n3 = simplify(C3-B3)/le;

% 连线矢量和平面夹角
sin_gama12 = dot(n1,C1C2)/d_C1C2;
sin_gama21 = dot(n2,-C1C2)/d_C1C2;
sin_gama23 = dot(n2,C2C3)/d_C2C3;
sin_gama32 = dot(n3,-C2C3)/d_C2C3;
sin_gama31 = dot(n3,C3C1)/d_C3C1;
sin_gama13 = dot(n1,-C3C1)/d_C3C1;

% 圆面线段投影
d_C1D12 = lm*sqrt(1-sin_gama12^2);
d_C2D21 = lm*sqrt(1-sin_gama21^2);
d_C2D23 = lm*sqrt(1-sin_gama23^2);
d_C3D32 = lm*sqrt(1-sin_gama32^2);
d_C3D31 = lm*sqrt(1-sin_gama31^2);
d_C1D13 = lm*sqrt(1-sin_gama13^2);

d12 = simplify( d_C1C2 - d_C1D12 - d_C2D21 );
d23 = simplify( d_C2C3 - d_C2D23 - d_C3D32 );
d31 = simplify( d_C3C1 - d_C3D31 - d_C1D13 );
% %test
% d12 = eval(d12);
% d23 = eval(d23);
% d31 = eval(d31);

%% SS管扭转角度delay
% 如按薄壁管计算 - 厚度<=半径/10
% J = pi*(ODs_i-IDs_i)*ODs_i^3/8;            %ss管极惯性矩, m^4

G = 77*10^9;                            % ss管剪切模量（奥氏体1Cr18Ni9Ti不锈钢），Pa
theta_max = deg2rad(3);                 % ss管最大允许扭转角度 5°

%tube 2
J_2 = pi*(ODs_2^4-IDs_2^4)/32;          % 极惯性矩，m^4
kz_2=G*J_2;                             % 扭转刚度，m^4*Pa
L_2 = delta_L_2+le+lb+lc;               % 管2总长度
theta_2 = tau_2*(L_2)/kz_2;               % 扭转角度,rad

%tube 1
IDs_1 = ODs_2+2*delta;
J_1 = pi*(ODs_1^4-IDs_1^4)/32;          % 极惯性矩， m^4
kz_1=G*J_1;                             % 扭转刚度，m^4*Pa
L_1 = delta_L_1+le+lb+lc;               % 管1总长度
theta_1 = tau_1*(L_1)/kz_1;               % 扭转角度,rad

% 如按非薄壁管计算
% J = pi*(ODs_i^4-IDs_i^4)/32;            %ss管极惯性矩, m^4


%% 优化问题 - 基于问题求解有约束非线性问题
%创建优化问题
prob = optimproblem;

%优化变量
ODs_1 = optimvar('ODs_1','LowerBound',ODn_1+2*b+2*w,'UpperBound',a);
ODs_2 = optimvar('ODs_2','LowerBound',max(ODn_2+2*b+2*w,ODn_1+2*b-2*delta),'UpperBound',a-2*w-2*delta);
le = optimvar('le','LowerBound',le_min,'UpperBound',1);
lb = optimvar('lb','LowerBound',0,'UpperBound',1);
alpha = optimvar('alpha','LowerBound',0,'UpperBound',pi/2);
beta_1 = optimvar('beta_1','LowerBound',0,'UpperBound',pi);
beta_2 = optimvar('beta_2','LowerBound',0,'UpperBound',pi);


%目标函数
%prob.Objective = -((11*lb)/50 + (11*le)/50 + 143/5000)/(2406250000*pi*((ODs_2 + 1/5000)^4 - ODs_1^4))  %薄壁内管扭转角最小
prob.Objective = le+lb    % 尺寸最compact
prob.ObjectiveSense = 'minimize';

%约束条件
cons_1 = ODs_2<= ODs_1-2*w-2*delta;             %SS管外径ODs_2的尺寸约束
cons_2 = (ODs_1*alpha)/(2*lb)<=epsilon_max;     %线性应变约束
% theta_1==theta_2
cons_3 = -((11*lb)/50 + (11*le)/50 + 143/5000)/(2406250000*pi*((ODs_2 + 1/5000)^4 - ODs_1^4))==(lb/5 + le/5 + 19/500)/(2406250000*pi*(ODs_2^4 - 4111825577611637/77371252455336267181195264)); %最内SS管扭转角约束
cons_4 = -((11*lb)/50 + (11*le)/50 + 143/5000)/(2406250000*pi*((ODs_2 + 1/5000)^4 - ODs_1^4))<=theta_max;
%collision avoidance
cons_5 = ((sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/2500)^2 + (le*sin(alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - (lb*(cos(alpha) - 1))/alpha + 7188907043224131/1152921504606846976)^2)^(1/2) - (3*(1 - (sin(alpha)^2*(le*sin(alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - (lb*(cos(alpha) - 1))/alpha + 7188907043224131/1152921504606846976)^2)/((sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/2500)^2 + (le*sin(alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - (lb*(cos(alpha) - 1))/alpha + 7188907043224131/1152921504606846976)^2))^(1/2))/250 - (3*(1 - (((cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + (lb*cos(beta_1)*(cos(alpha) - 1))/alpha)*(le*sin(alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - (lb*(cos(alpha) - 1))/alpha + 7188907043224131/1152921504606846976))/le - ((sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/2500)*(sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + (lb*sin(beta_1)*(cos(alpha) - 1))/alpha))/le)^2/((sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/2500)^2 + (le*sin(alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - (lb*(cos(alpha) - 1))/alpha + 7188907043224131/1152921504606846976)^2))^(1/2))/250>=d_min;
cons_6 = ((cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha))^2 + (sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/1250)^2)^(1/2) - (3*(1 - (cos(beta_1 + beta_2)*sin(alpha)*(cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha)) - sin(beta_1 + beta_2)*sin(alpha)*(sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/1250))^2/((cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha))^2 + (sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/1250)^2))^(1/2))/250 - (3*(1 - (cos(beta_1)*sin(alpha)*(cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha)) - sin(alpha)*sin(beta_1)*(sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/1250))^2/((cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - cos(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha))^2 + (sin(beta_1)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) + 9/1250)^2))^(1/2))/250>=d_min;
cons_7 = ((sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - 9/2500)^2 + (cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - le*sin(alpha) + (lb*(cos(alpha) - 1))/alpha - 7188907043224131/1152921504606846976)^2)^(1/2) - (3*(1 - (sin(alpha)^2*(cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - le*sin(alpha) + (lb*(cos(alpha) - 1))/alpha - 7188907043224131/1152921504606846976)^2)/((sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - 9/2500)^2 + (cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - le*sin(alpha) + (lb*(cos(alpha) - 1))/alpha - 7188907043224131/1152921504606846976)^2))^(1/2))/250 - (3*(1 - (sin(beta_1 + beta_2)*sin(alpha)*(sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - 9/2500) + cos(beta_1 + beta_2)*sin(alpha)*(cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - le*sin(alpha) + (lb*(cos(alpha) - 1))/alpha - 7188907043224131/1152921504606846976))^2/((sin(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - 9/2500)^2 + (cos(beta_1 + beta_2)*(le*sin(alpha) - (lb*(cos(alpha) - 1))/alpha) - le*sin(alpha) + (lb*(cos(alpha) - 1))/alpha - 7188907043224131/1152921504606846976)^2))^(1/2))/250>=d_min;

%在问题中包含约束 - 其余的约束直接在定义优化变量的时候在上下界处给出了
prob.Constraints.cons1 = cons_1;
prob.Constraints.cons2 = cons_2;
prob.Constraints.cons3 = cons_3;
prob.Constraints.cons4 = cons_4;
prob.Constraints.cons5 = cons_5;
prob.Constraints.cons3 = cons_6;
prob.Constraints.cons4 = cons_7;

%检查此优化问题
show(prob)

%初值
initialpt.ODs_1 = 0;
initialpt.ODs_2 = 0;
initialpt.le = 0;
initialpt.lb = 0;
initialpt.alpha = 0;
initialpt.beta_1 = 0;
initialpt.beta_2 = 0;

%求解问题
[sol,fval] = solve(prob,initialpt)

%检查解
nx = norm([sol.ODs_1 sol.ODs_2 sol.lb sol.le sol.alpha sol.beta_1 sol.beta_2])                        %范数要确保小于或等于1


%% 优化结果手动输出
alpha_deg = rad2deg(sol.alpha) 
theta_1 = -((11*sol.lb)/50 + (11*sol.le)/50 + 143/5000)/(2406250000*pi*((sol.ODs_2 + 1/5000)^4 - sol.ODs_1^4))
theta_2 = (sol.lb/5 + sol.le/5 + 19/500)/(2406250000*pi*(sol.ODs_2^4 - 4111825577611637/77371252455336267181195264))
theta_1_deg = rad2deg(theta_1)
theta_2_deg = rad2deg(theta_2)
 
IDs_2
ODs_2 = sol.ODs_2
IDs_1 = ODs_2+2*delta
ODs_1 = sol.ODs_1



