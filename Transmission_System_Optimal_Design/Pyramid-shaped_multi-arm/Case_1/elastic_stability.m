%% Elastic Stability of "SS + NiTi" tube pair
% by Chao Zhang
% Date：2022/10/20

clear;clc;


%% Load CTR Params
CTR_params;


%% NiTi tube params
% Outer tube
k1x = pi*En*(ODn_1^4-IDn_1^4)/64;       % Bending stiffness
k1y = k1x;
k1z = k1x/(1+poisson_rate);             % Torsional stiffness
hat_u1 = [1/r1;0;0];                    % Pre-curvature,m^(-1)

% Inner tube
k2x = pi*En*(ODn_2^4-IDn_2^4)/64;
k2y = k2x;
k2z = k2x/(1+poisson_rate);
hat_u2 = [1/r2;0;0];


%% Design Results of SS tubes
% Empirically designed result
eODs_1 = 4.0e-3;
eIDs_1 = 3.3e-3;
eODs_2 = 3.0e-3;
eIDs_2 = 2.4e-3;
eL_1 = 500e-3;
eL_2 = 545e-3;
ek1z_s = Gs*pi*(eODs_1^4-eIDs_1^4)/32;       % Torsional stiffness of SS tube 1
ek2z_s = Gs*pi*(eODs_2^4-eIDs_2^4)/32;       % Torsional stiffness of SS tube 2

% Optimally designed result of point A
oaODs_1 = 5e-3;
oaIDs_1 = 4.39e-3;
oaODs_2 = 4.19e-3;
oaIDs_2 = 2.4e-3;
oaL_1 = 524.8e-3;
oaL_2 = 569.8e-3;
oak1z_s = Gs*pi*(oaODs_1^4-oaIDs_1^4)/32;
oak2z_s = Gs*pi*(oaODs_2^4-oaIDs_2^4)/32;

% Optimally designed result of point B
obODs_1 = 4.05e-3;
obIDs_1 = 3.63e-3;
obODs_2 = 3.43e-3;
obIDs_2 = 2.4e-3;
obL_1 = 471.8e-3;
obL_2 = 516.8e-3;
obk1z_s = Gs*pi*(obODs_1^4-obIDs_1^4)/32;
obk2z_s = Gs*pi*(obODs_2^4-obIDs_2^4)/32;

% Optimally designed result of point C
ocODs_1 = 4.8e-3;
ocIDs_1 = 4.2e-3; 
ocODs_2 = 4.0e-3;
ocIDs_2 = 2.4e-3;
ocL_1 = 512.8e-3;
ocL_2 = 557.8e-3;
ock1z_s = Gs*pi*(ocODs_1^4-ocIDs_1^4)/32;
ock2z_s = Gs*pi*(ocODs_2^4-ocIDs_2^4)/32;


%% Computing the Elastic Stability
% 弹性作用扭矩
c = (1+poisson_rate)*norm(hat_u1)*norm(hat_u2);     % Constant c in Dupont 2010 T-RO eq.(24)
L_sqrtc = s1*sqrt(c);                               % Constant L*sqrt(c) in Dupont 2010 T-RO eq.(34)

alpha_L = linspace(0,2*pi,200);                     % Relative rotation angle of tube pair at distal end
alpha_0 = zeros(length(alpha_L),1);                 % Relatice rotation angle of tube pair at proximal end
tau = zeros(length(alpha_L),1);                     % Generated torque
tors_comp = zeros(length(alpha_L),12); 

for i=1:length(alpha_L)
    % Jacobi elliptic function in Dupont 2010 T-RO eq.(34)
    alpha_0(i) = 2*acos(cos(alpha_L(i)/2)*jacobiCD(L_sqrtc,cos(alpha_L(i)/2)^2));
    
    % Generated torque at the base of the NiTi tube pair
    if (alpha_L(i) <= pi)
        tau(i) = k1z*k2z/(k1z+k2z)*sqrt(2*c*(cos(alpha_L(i))-cos(alpha_0(i))));
    else
        tau(i) = -k1z*k2z/(k1z+k2z)*sqrt(2*c*(cos(alpha_L(i))-cos(alpha_0(i))));
    end
    
    % Calculate the relative twist angle
    ezeta_1 = tau(i)*eL_1/ek1z_s;           % Empirically design - twist angle of SS tube 1
    ezeta_2 = -tau(i)*eL_2/ek2z_s;          % Empirically design - twist angle of SS tube 2
    oazeta_1 = tau(i)*oaL_1/oak1z_s;        % Optimally design (A) - twist angle of SS tube 1
    oazeta_2 = -tau(i)*oaL_2/oak2z_s;       % Optimally design (A) - twist angle of SS tube 2
    obzeta_1 = tau(i)*obL_1/obk1z_s;        % Optimally design (B) - twist angle of SS tube 1
    obzeta_2 = -tau(i)*obL_2/obk2z_s;       % Optimally design (B) - twist angle of SS tube 2
    oczeta_1 = tau(i)*ocL_1/ock1z_s;        % Optimally design (C) - twist angle of SS tube 1
    oczeta_2 = -tau(i)*ocL_2/ock2z_s;       % Optimally design (C) - twist angle of SS tube 2
    nzeta_1 = tau(i)*obL_1/k1z;             % Only NiTi tube - twist angle of SS tube 1
    nzeta_2 = -tau(i)*obL_2/k2z;            % Only NiTi tube - twist angle of SS tube 2
    
    % Record the values
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
tau_1_max = max(tau_1)             % Maximum generated torque in outer tube
tau_2_max = max(tau_2)             % Maximum generated torque in inner tube
ezeta_1 = max(tors_comp(:,3))
ezeta_2 = max(tors_comp(:,4))


%% Plot the Rotation angles at the Poximal and Distal end of Tube 2 (alpha_ss(0) vs. alpha_L)
figure('Name','alpha_ss vs alpha_L');

% Using totally rigid SS tube
plot(tors_comp(:,2),tors_comp(:,1),'--','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880]);
hold on

% Optimally designed SS tube - point A
plot(tors_comp(:,2)+tors_comp(:,5)-tors_comp(:,6),tors_comp(:,1),'-','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560]);
hold on

% Optimally designed SS tube - point B
plot(tors_comp(:,2)+tors_comp(:,7)-tors_comp(:,8),tors_comp(:,1),'-b','LineWidth',1.5);
hold on

% Optimally designed SS tube - point C
plot(tors_comp(:,2)+tors_comp(:,9)-tors_comp(:,10),tors_comp(:,1),'-c','LineWidth',1.5);
hold on

% Empirically designed SS tube
plot(tors_comp(:,2)+tors_comp(:,3)-tors_comp(:,4),tors_comp(:,1),'-r','LineWidth',1.5);
hold on

% Only NiTi tube
%plot(tors_comp(:,2)+tors_comp(:,11)-tors_comp(:,12),tors_comp(:,1),'--m','LineWidth',1.5);
%hold on

line([180 180],[0 360],'Color','black','LineStyle','-.');
axis equal;
grid on;
set(gca,'XTick',0:30:360);
set(gca,'YTick',0:30:360);
xlabel('\alpha_S(0) (\circ)');
ylabel('\alpha_N(s_1) (\circ)');
% Unit \degree
xlim([0 360]);
ylim([0 360]);
legend('Rigid SS tube','Optimal design at point A','Optimal design at point B','Optimal design at point C','Empirical design');

% % 导出设置
% %大小 - 20x20cm
% %渲染 - 600dpi
% %字体 - 14磅，Times New Roman



