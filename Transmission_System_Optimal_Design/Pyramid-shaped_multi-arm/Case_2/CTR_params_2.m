%% Design Parameters of the NiTi Manipulator
% It is a four-armed CTR with 2 tubes in each arm
% by Chao Zhang
% Date：2022/10/20


%% Parameters of NiTi Manipulator
n = 4;                                  % Number of arms
a = 5e-3;                               % Center distrance between adjacent arms, m
D = 12e-3;                              % Straight sheath diameter

En = 51.9e9;                            % NiTi elastic modulus, Pa
poisson_rate = 0.3;                     % NiTi Poisson's ratio
Gn = En/(2*(1+poisson_rate));           % NiTi shear modulus

ODn_1 = 1.8e-3;                         % Outer diameter of NiTi tube 1 (outer tube)
IDn_1 = 1.46e-3;                        % Inner diameter of NiTi tube 1 (outer tube)
ODn_2 = 1.26e-3;                        % Outer diameter of NiTi tube 2 (inner tube)
IDn_2 = 1.1e-3;                         % Inner diameter of NiTi tube 2 (inner tube)

r1 = 120e-3;                            % Pre-curvature of tube 1
s1 = 60e-3;                             % Arc lenth of tube 1
r2 = r1;                                % Pre-curvature of tube 2
s2 = s1;                                % Arc lenth of tube 2

tau_1_max = 0.00951263417855761;        % Generated torque, calculated by "elastic_stability_2.m"
tau_2_max = 0.00951263417855761;


%% Parameters of Stainless Steel (SS) Tubes (SUS 316L)
Es = 193e9;                             % SS elastic modulus, Pa 
Gs = 72e9;                              % SS shear modulus
yeild_s = 310e6;                        % SS tensile yield stress
epsilon_e = yeild_s/Es;                 % SS elastic strain limit

zeta_max = deg2rad(3);                  % SS maximum allowable relatice twist angle, deg

b = 0.1e-3;                             % Minimum bonding thickness bewtween SS and NiTi tubes
w = 0.1e-3;                             % Minimum wall thickness of SS tube
delta = 0.1e-3;                         % Minimum clearence between SS tubes

lc = 80e-3;                             % Length of the straight sheath
le_min = 40e-3;                         % Minimum allowable le
delta_L_1 = 60e-3;                  	% Distance from the first actuation module forward limiting position to the bending unit
delta_L_2 = 115e-3;                     % Distance from the second actuation module forward limiting position to the bending unit
ra = 18e-3;                             % Radius of collision hazard region