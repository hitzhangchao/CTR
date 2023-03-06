%% 该对NiTi balanced pair的弹性稳定性分析
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

% 弹性作用扭矩
c = (1+poisson_rate)*norm(hat_u1)*norm(hat_u2);     % Dupont公式(24)中的常数c
L_sqrtc = s1*sqrt(c);                               % Dupont公式(34)中的常数L*sqrt(c)

alpha_L = linspace(0,2*pi,200);                     % balance pair末端相对扭转角
alpha_0 = zeros(length(alpha_L),1);                 % balance pair近端相对扭转角
tau = zeros(length(alpha_L),1);                     % 弹性作用扭矩
for i=1:length(alpha_L)
    alpha_0(i) = 2*acos(cos(alpha_L(i)/2)*jacobiCD(L_sqrtc,cos(alpha_L(i)/2)^2));% Jacobi椭圆函数，公式(34)
    
    % 计算NiTi管近端扭矩(因为有一个开根号，正负需定号)
    if (alpha_L(i) <= pi)
        tau(i) = k1z*k2z/(k1z+k2z)*sqrt(2*c*(cos(alpha_L(i))-cos(alpha_0(i))));
    else
        tau(i) = -k1z*k2z/(k1z+k2z)*sqrt(2*c*(cos(alpha_L(i))-cos(alpha_0(i))));
    end
end
tau_1 = tau;
tau_2 = -tau;
tau_1_max = max(tau_1)             % 外管最大扭矩
tau_2_max = max(tau_2)             % 内管最大扭矩
% 提前计算一下最大相对扭转角
% L1 = 450e-3;                % 管1总长度
% L2 = 500e-3;                % 管2总长度
% ODs_1 = 4e-3;
% IDs_1 = 3.6e-3;
% ODs_2 = 3.2e-3;
% IDs_2 = 2.7e-3;
% zeta_1_max = rad2deg(32*tau_1_max*L1/(pi*Gs*(ODs_1^4-IDs_1^4)))
% zeta_2_max = rad2deg(32*tau_2_max*L2/(pi*Gs*(ODs_2^4-IDs_2^4)))

% 绘制NiTi管alpha_0 和 alpha_L 关系
figure('Name','alpha_0 vs alpha_L');
plot(rad2deg(alpha_0),rad2deg(alpha_L),'-r','LineWidth',2);
hold on
line([180 180],[0 360],'Color','black','LineStyle','-.');
set(gca,'XTick',[0:30:360]);
set(gca,'YTick',[0:30:360]);
xlim([0 360]);
ylim([0 360]);
axis equal;
grid on;
xlabel('\alpha_0 [rad]');
ylabel('\alpha_L [rad]');

% 绘制NiTi管alpha_0 和 tau 关系
figure('Name','alpha_0 vs tau');
plot(rad2deg(alpha_0),tau_1,'-b','LineWidth',2);
%hold on
%plot(rad2deg(alpha_0),tau_2,'-m','LineWidth',2);
set(gca,'XTick',[0:30:360]);
xlim([0 360]);
grid on;
xlabel('\alpha_0 [rad]');
ylabel('\tau [Nm]');

