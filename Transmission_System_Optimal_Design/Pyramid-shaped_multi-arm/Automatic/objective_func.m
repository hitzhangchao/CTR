function F = objective_func(x)
%x(1);      %ODs_1
%x(2);      %ODs_2
%x(3);      %lb
%x(4);      %alpha
%x(5);      %le

%x(6);      %beta_1
%x(7);      %beta_2

%% Load CTR params
CTR_params;

%% Geometric Constraints - Relations of tube diameters
IDs_2 = ODn_2+2*b;                      % 最内SS管内径

%% Geometric Constraints - Shape and dimensions of SS tubes
L2 = delta_L_2+x(5)+x(3)+lc;                % 管2总长度

%% Deformation Constraints - Torsional twisting of SS tubes
% 内管
J2 = pi*(x(2)^4-IDs_2^4)/32;           % 极惯性矩，m^4
kz_2 = Gs*J2;                           % 扭转刚度，m^4*Pa
zeta_2 = tau_2_max*L2/kz_2;           % 扭转角度,rad

%% function
F1 = (180/pi)*zeta_2;               % SS管相对扭转角度，deg
F2 = x(3)+x(5);                     % SS管长

F = [F1,F2];
end
