%% 该对NiTi balanced pair + SS tube的弹性稳定性分析
% This case demonstrates the design optimization of a triple-arm CTR with 2 tubes (balanced pair) in each arm.
% by Zhang Chao
% Date：2022/10/20

clear;clc;

%% 加载NiTi管基础物理参数
CTR_params;

%% NiTi参数运算
% 外管
k1x = pi*En*(ODn_1^4-IDn_1^4)/64;        % 弯曲刚度
k1y = k1x;
k1z = k1x/(1+poisson_rate);             % 扭转刚度
hat_u1 = [1/r1;0;0];                    % 预曲率,m^(-1)

% 内管
k2x = pi*En*(ODn_2^4-IDn_2^4)/64;
k2y = k2x;
k2z = k2x/(1+poisson_rate);
hat_u2 = [1/r2;0;0];


%% SS管设计结果
% 经验设计结果
eODs_1 = 4.0e-3;
eIDs_1 = 3.3e-3;
eODs_2 = 3.0e-3;
eIDs_2 = 2.4e-3;
eL_1 = 500e-3;
eL_2 = 545e-3;
ek1z_s = Gs*pi*(eODs_1^4-eIDs_1^4)/32;       % SS管1的扭转刚度
ek2z_s = Gs*pi*(eODs_2^4-eIDs_2^4)/32;       % SS管2的扭转刚度

% 优化设计的结果 （point A）
oaODs_1 = 5e-3;
oaIDs_1 = 4.39e-3;
oaODs_2 = 4.19e-3;
oaIDs_2 = 2.4e-3;
oaL_1 = 524.8e-3;
oaL_2 = 569.8e-3;
oak1z_s = Gs*pi*(oaODs_1^4-oaIDs_1^4)/32;       % SS管1的扭转刚度
oak2z_s = Gs*pi*(oaODs_2^4-oaIDs_2^4)/32;       % SS管2的扭转刚度

% 优化设计的结果 （point B）
obODs_1 = 4.05e-3;
obIDs_1 = 3.63e-3;
obODs_2 = 3.43e-3;
obIDs_2 = 2.4e-3;
obL_1 = 471.8e-3;
obL_2 = 516.8e-3;
obk1z_s = Gs*pi*(obODs_1^4-obIDs_1^4)/32;       % SS管1的扭转刚度
obk2z_s = Gs*pi*(obODs_2^4-obIDs_2^4)/32;       % SS管2的扭转刚度

% 优化设计的结果 （point C）
ocODs_1 = 4.8e-3;
ocIDs_1 = 4.2e-3; 
ocODs_2 = 4.0e-3;
ocIDs_2 = 2.4e-3;
ocL_1 = 512.8e-3;
ocL_2 = 557.8e-3;
ock1z_s = Gs*pi*(ocODs_1^4-ocIDs_1^4)/32;       % SS管1的扭转刚度
ock2z_s = Gs*pi*(ocODs_2^4-ocIDs_2^4)/32;       % SS管2的扭转刚度


%% 计算弹性稳定性
% 弹性作用扭矩
c = (1+poisson_rate)*norm(hat_u1)*norm(hat_u2);     % Dupont公式(24)中的常数c
L_sqrtc = s1*sqrt(c);                               % Dupont公式(34)中的常数L*sqrt(c)

alpha_L = linspace(0,2*pi,200);                     % balance pair末端相对扭转角
alpha_0 = zeros(length(alpha_L),1);                 % balance pair近端相对扭转角
tau = zeros(length(alpha_L),1);                     % 弹性作用扭矩
tors_comp = zeros(length(alpha_L),12);              % 预分配存储计算结果的数组空间

for i=1:length(alpha_L)
    % Jacobi椭圆函数，Dupont论文公式(34)
    alpha_0(i) = 2*acos(cos(alpha_L(i)/2)*jacobiCD(L_sqrtc,cos(alpha_L(i)/2)^2));
    
    % 计算NiTi管近端扭矩(因为有一个开根号，正负需定号)
    if (alpha_L(i) <= pi)
        tau(i) = k1z*k2z/(k1z+k2z)*sqrt(2*c*(cos(alpha_L(i))-cos(alpha_0(i))));
    else
        tau(i) = -k1z*k2z/(k1z+k2z)*sqrt(2*c*(cos(alpha_L(i))-cos(alpha_0(i))));
    end
    
    % 计算扭转角等
    ezeta_1 = tau(i)*eL_1/ek1z_s;           % 经验设计 - SS管1扭转角delay
    ezeta_2 = -tau(i)*eL_2/ek2z_s;          % 经验设计 - SS管2扭转角delay
    oazeta_1 = tau(i)*oaL_1/oak1z_s;        % 优化设计（A点） - SS管1扭转角delay
    oazeta_2 = -tau(i)*oaL_2/oak2z_s;       % 优化设计（A点） - SS管2扭转角delay
    obzeta_1 = tau(i)*obL_1/obk1z_s;        % 优化设计（B点） - SS管1扭转角delay
    obzeta_2 = -tau(i)*obL_2/obk2z_s;       % 优化设计（B点） - SS管2扭转角delay
    oczeta_1 = tau(i)*ocL_1/ock1z_s;        % 优化设计（C点） - SS管1扭转角delay
    oczeta_2 = -tau(i)*ocL_2/ock2z_s;       % 优化设计（C点） - SS管2扭转角delay
    nzeta_1 = tau(i)*obL_1/k1z;             % 纯NiTi管 - SS管1扭转角delay
    nzeta_2 = -tau(i)*obL_2/k2z;            % 纯NiTi管 - SS管2扭转角delay
    
    %单位用 degree
    tors_comp(i,1) = rad2deg(alpha_L(i));
    tors_comp(i,2) = rad2deg(alpha_0(i));
    tors_comp(i,3) = rad2deg(ezeta_1);
    tors_comp(i,4) = rad2deg(ezeta_2);
    tors_comp(i,5) = rad2deg(oazeta_1);
    tors_comp(i,6) = rad2deg(oazeta_2);
    tors_comp(i,7) = rad2deg(obzeta_1);
    tors_comp(i,8) = rad2deg(obzeta_2);
    tors_comp(i,9) = rad2deg(oczeta_1);
    tors_comp(i,10) = rad2deg(oczeta_2);
    tors_comp(i,11) = rad2deg(nzeta_1);
    tors_comp(i,12) = rad2deg(nzeta_2);
end
tau_1 = tau;
tau_2 = -tau;
tau_1_max = max(tau_1)             % 外管最大扭矩
tau_2_max = max(tau_2)             % 内管最大扭矩
ezeta_1 = max(tors_comp(:,3))
ezeta_2 = max(tors_comp(:,4))


%% 绘制SS管2近端旋转和NiTi管2远端转动关系 alpha_ss alpha_L 关系
figure('Name','alpha_ss vs alpha_L');
% 假设SS管2完全刚性
plot(tors_comp(:,2),tors_comp(:,1),'--','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880]);
hold on
% 优化设计 - point A
plot(tors_comp(:,2)+tors_comp(:,5)-tors_comp(:,6),tors_comp(:,1),'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560]);
hold on
% 优化设计 - point B
plot(tors_comp(:,2)+tors_comp(:,7)-tors_comp(:,8),tors_comp(:,1),'-b','LineWidth',1.5);
hold on
% 优化设计 - point C
plot(tors_comp(:,2)+tors_comp(:,9)-tors_comp(:,10),tors_comp(:,1),'-c','LineWidth',1.5);
hold on
% 经验设计
plot(tors_comp(:,2)+tors_comp(:,3)-tors_comp(:,4),tors_comp(:,1),'-r','LineWidth',1.5);
hold on
% 假设纯采用NiTi管
%plot(tors_comp(:,2)+tors_comp(:,11)-tors_comp(:,12),tors_comp(:,1),'--m','LineWidth',1.5);
%hold on

line([180 180],[0 360],'Color','black','LineStyle','-.');
axis equal;
grid on;
set(gca,'XTick',[0:30:360]);
set(gca,'YTick',[0:30:360]);
xlabel('\alpha_S(0) (\circ)');
ylabel('\alpha_N(s_1) (\circ)');
% 单位用degree
xlim([0 360]);
ylim([0 360]);
legend('Rigid SS tube','Optimal design at point A','Optimal design at point B','Optimal design at point C','Empirical design');
% % 导出设置
% %大小 - 20x20cm
% %渲染 - 600dpi
% %字体 - 14磅，Times New Roman



