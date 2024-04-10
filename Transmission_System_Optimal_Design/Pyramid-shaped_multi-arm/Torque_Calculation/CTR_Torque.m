%% Maximum Generated Torque Calculation - Three NiTi Tubes
% by Zhang Chao
% Date：2023/9/28

clear;clc;

%% Load Parameters
CTR_params;


%% Sampling set
sample_cnt = 1000;                  % Numbers of Sampling
sample_Rp = zeros(6+3,sample_cnt);  % Matrix to storage 6 joint variablese and 3 torque values

rng('shuffle');
tic
for sp_cnt=1:sample_cnt
    % Kinematic Input
    theta_1 = 0;                % Tube 1 (outermost) does not rotate
    theta_2 = -pi+2*pi*rand;    % Rotation angle of tube 2 (middle)
    theta_3 = -pi+2*pi*rand;    % Rotation angle of tube 3 (innermost)
    len_ex_1 = 30e-3;           % Tube 1 fully extended
    len_ex_2 = 30.01e-3;        % Tube 2 fully retracted
    len_ex_3 = 30.02e-3;        % Tube 2 fully retracted

    theta = [theta_1,theta_2,theta_3];
    len_ex = [len_ex_1,len_ex_2,len_ex_3];

    % Forward Kinematics Computing using Torsionally Compliant Model
    [p,T_tip,d_tip,ss,theta_dot_0] = ctr_fk_compliant(len_ex,theta,len,len_cu,uy_star,I);

    % Calculate the torque in the base of eath NiTi tube
    kz = G*J';                  % Tube stiffness along the central axis
    tau = kz .* theta_dot_0;    % Generated torque
    
    % Data storage
    sample_Rp(1:6,sp_cnt) = [theta,len_ex]';
    sample_Rp(7:9,sp_cnt) = tau;
end
toc

%% Maximum Generated Torque
% Tube 1
tau_1 = sample_Rp(7,:);
max_tau_1 = max(tau_1)
index_tau_1 = find(tau_1==max(tau_1));

% Tube 2
tau_2 = sample_Rp(8,:);
max_tau_2 = max(tau_2)
index_tau_2 = find(tau_2==max(tau_2));

% Tube 3
tau_3 = sample_Rp(9,:);
max_tau_3 = max(tau_3)
index_tau_3 = find(tau_3==max(tau_3));



