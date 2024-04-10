%% Design Parameters of the NiTi Manipulator
% It is a Four-armed CTR with 3 tubes in each arm
% by Chao Zhang
% Date：2023/10/1


%% Parameters about Tubes
% NiTi Tubes
En = 51.9e9;                            % NiTi elastic modulus, Pa
poisson_rate = 0.3;                     % NiTi Poisson's ratio
Gn = En/(2*(1+poisson_rate));           % NiTi shear modulus

ODn_1 = 2.25e-3;                        % Outer diameter of NiTi tube 1 (outermost tube), m
IDn_1 = 1.75e-3;                        % Inner diameter of NiTi tube 1 (outermost tube)
ODn_2 = 1.65e-3;                        % Outer diameter of NiTi tube 2 (middle tube)
IDn_2 = 1.40e-3;                        % Inner diameter of NiTi tube 2 (middle tube)
ODn_3 = 1.26e-3;                        % Outer diameter of NiTi tube 3 (innermost tube)
IDn_3 = 1.1e-3;                         % Inner diameter of NiTi tube 3 (innermost tube)

r1 = 135e-3;                            % Pre-curvature of tube 1, m
s1 = 45e-3;                             % Arc lenth of tube 1
r2 = 85e-3;
s2 = 80e-3;
r3 = 65e-3;
s3 = 110e-3;

% SS Tubes
Es = 193e9;                             % SS elastic modulus, Pa
Gs = 72e9;                              % SS shear modulus (SUS 316L)
yeild_s = 310e6;                        % SS tensile yield stress
epsilon_e = yeild_s/Es;                 % SS elastic strain limit
zeta_max = deg2rad(4.5);                % SS maximum allowable relatice twist angle, deg

% Bonding parameters of NiTi and SS tubes
b = 0.1e-3;                             % Minimum bonding thickness bewtween SS and NiTi tubes
w = 0.1e-3;                             % Minimum wall thickness of SS tube
delta = 0.1e-3;                         % Minimum clearence between SS tubes


%% Parameters about Whole System
n = 4;                                  % Number of arms
a = 5e-3;                               % Center distance between adjacent arms, m
D = 12e-3;                              % Straight sheath diameter
ls = 80e-3;                             % Length of straight sheath
la = 15e-3;                             % Minimum assembly distance between the distal ends of adjacent SS tubes
le_min = 40e-3;                         % Minimum allowable le
delta_L_1 = 100e-3;                  	% Distance from the first actuation module forward limiting position to the bending unit
delta_L_2 = 145e-3;                     % Distance from the second actuation module forward limiting position to the bending unit
delta_L_3 = 170e-3;                     % Distance from the third actuation module forward limiting position to the bending unit
ra = 18e-3;                             % Radius of collision hazard region

tau_1_max = 0.0247;          % Generated torque, calculated by ".../Torque_Calculation/CTR_Torque.m"
tau_2_max = 0.0202;
tau_3_max = 0.0093;
