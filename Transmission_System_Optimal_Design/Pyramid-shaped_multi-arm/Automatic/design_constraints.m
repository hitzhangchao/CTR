function [C,Ceq] = design_constraints(x)
%x(1);      %ODs_1
%x(2);      %ODs_2
%x(3);      %lb
%x(4);      %alpha
%x(5);      %le

%x(6);      %beta_1
%x(7);      %beta_2

syms beta_1 real

%% Load CTR params
CTR_params;

%% Geometric Constraints - Relations of tube diameters
IDs_2 = ODn_2+2*b;                      % 最内SS管内径

%% Geometric Constraints - Shape and dimensions of SS tubes
R = x(3)/x(4);                           % 圆弧段弯曲半径
L1 = delta_L_1+x(5)+x(3)+lc;                % 管1总长度
L2 = delta_L_2+x(5)+x(3)+lc;                % 管2总长度

%% Geometric Constraints - Collision avoidance of actuation units
% 建立坐标系，空间各点坐标
B1 = [sqrt(3)/3*a+R*(1-cos(x(4)));0;lc+R*sin(x(4))];
B2 = [-sqrt(3)/6*a+R*(1-cos(x(4)))*cos(beta_1);a/2+R*(1-cos(x(4)))*sin(beta_1);lc+R*sin(x(4))];

C = [0;0;lc+R*sin(x(4))+x(5)*cos(x(4))];
C1 = [sqrt(3)/3*a+R*(1-cos(x(4)))+x(5)*sin(x(4));0;lc+R*sin(x(4))+x(5)*cos(x(4))];
C2 = [-sqrt(3)/6*a+(R*(1-cos(x(4)))+x(5)*sin(x(4)))*cos(beta_1);a/2+(R*(1-cos(x(4)))+x(5)*sin(x(4)))*sin(beta_1);lc+R*sin(x(4))+x(5)*cos(x(4))];

% actuation module中心点连线矢量
C1C2 = C2-C1;
d_C1C2 = sqrt(C1C2(1)^2+C1C2(2)^2);

% actuation module平面法向量
n1 = (C1-B1)/x(5);
n2 = (C2-B2)/x(5);

% 连线矢量和平面夹角
sin_gama12 = dot(n1,C1C2)/d_C1C2;
sin_gama21 = dot(n2,-C1C2)/d_C1C2;

% 圆面线段投影
d_C1D12 = ra*sqrt(1-sin_gama12^2);
d_C2D21 = ra*sqrt(1-sin_gama21^2);

d1 = d_C1C2 - d_C1D12 - d_C2D21;


%% Deformation Constraints - Elastic strain limit of bending
r = x(1)/2;                            % 最外SS管半径
epsilon = r/R;                          % SS管应变
%epsilon <= epsilon_e                    %线性应变约束


%% Deformation Constraints - Torsional twisting of SS tubes
% 内管
J2 = pi*(x(2)^4-IDs_2^4)/32;           % 极惯性矩，m^4
kz_2 = Gs*J2;                           % 扭转刚度，m^4*Pa
zeta_2 = tau_2_max*(L2)/kz_2;           % 扭转角度,rad

% 外管
IDs_1 = x(2)+2*delta;                  % 外管内径
J1 = pi*(x(1)^4-IDs_1^4)/32;           % 极惯性矩， m^4
kz_1 = Gs*J1;                           % 扭转刚度，m^4*Pa
zeta_1 = tau_1_max*(L1)/kz_1;           % 扭转角度,rad

%% design constraints
% tube diameter -- ODs_2<= ODs_1-2*w-2*delta
C(1) = x(2) - x(1) + 2*w + 2*delta;

%epsilon <= epsilon_e
C(2) = epsilon - epsilon_e;

%zeta_i<=zeta_max，单位deg
C(3) = (180/pi)*(zeta_1 - zeta_max);

%d1>=0 -> -d1<=0  (我们令beta=2*pi/n)
C(4) = - subs(d1,beta_1,2*pi/n);

%zeta_1=zeta_2
Ceq(1) = (180/pi)*(zeta_1 - zeta_2);
end

