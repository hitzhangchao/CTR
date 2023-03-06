%% 三臂CTR，每臂由一对 NiTi balanced pair组成  - 设计参数
% 只是 NiTi manipulator
% by Zhang Chao
% Date：2022/10/20

%% 管子相关
% NiTi管相关
En = 48.7e9;                              % NiTinol材料弹性模量，Pa
poisson_rate = 0.3;                     % Nitinol材料泊松比
Gn = En/(2*(1+poisson_rate));           % Nitinol材料剪切模量,Pa

ODn_1 = 2.7e-3;                         % 外NiTi管外径
IDn_1 = 2.4e-3;                         % 外NiTi管内径
ODn_2 = 2.2e-3;                         % 内NiTi管外径
IDn_2 = 1.4e-3;                         % 内NiTi管内径

r1 = 135e-3;                            % 外NiTi管预弯曲半径
s1 = 160e-3;                            % 外NiTi管弧长
r2 = r1;                                % 外NiTi管预弯曲半径
s2 = s1;                                % 外NiTi管弧长

% SS管相关
Es = 193e9;                             % SS材料弹性模量，Pa
Gs = 72e9;                              % SS材料剪切模量（奥氏体1Cr18Ni9Ti不锈钢）
yeild_s = 310e6;                        % SS材料tensile yield stress（屈服应力）
epsilon_e = yeild_s/Es;                 % SS管线性应变段最大允许应变（如果超过，管子就可能塑性变形）
zeta_max = deg2rad(8);                  % SS管最大允许扭转角度

% NiTi管和SS粘接参数
b = 0.1e-3;                             % SS管和NiTi管最小粘接间隙
w = 0.1e-3;                             % SS管最小壁厚
delta = 0.1e-3;                         % SS管间最小游动间隙


%% 整机相关
n = 3;                                  % 手臂数
a = 5e-3;                             % 臂间中心距,m
D = 12e-3;                              % 单孔手术的孔径
lc = 80e-3;                             % SS即管伸出直线段的长度,准直器
le_min = 40e-3;                         % le最小距离
delta_L_1 = 100e-3;                  	% SS管1电机端到lp起点的距离
delta_L_2 = 145e-3;                     % SS管2电机端到le起点的距离
ra = 18e-3;                             % 驱动端安全圆半径（最大直径为35mm）

% s = 160e-3时
tau_1_max = 0.200283940448042;          % NiTi tube pair弹性作用扭矩，根据elastic_stability.m
tau_2_max = 0.200283940448042;