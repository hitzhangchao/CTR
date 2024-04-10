%% Design Parameters of the NiTi Manipulator
% by Zhang Chao
% Date：2022/9/23

n = 3;                      % number of tubes
E = 60E9;                   % NiTi elastic modulus, Pa
poisson_rate = 0.3;         % NiTi Poisson's ratio
G = E/(2*(1+poisson_rate)); % NiTi shear modulus


%% NiTi Tube 1 - Outermost
OD_1 = 2.25e-3;             % Outer diameter of NiTi tube 1, m
ID_1 = 1.75e-3;             % Inner diameter of NiTi tube 1 
len_1 = 40e-3;              % Total length of NiTi tube 1
len_cu_1 = 30e-3;           % Curved length of NiTi tube 1
r1 = 135e-3;                % Radius of curvature
u_y1_star = 1/r1;           % Pre-curvature (y-direction)
I1 = pi*(OD_1^4-ID_1^4)/64; % Cross-sectional moment of inertia, m^4
J1 = pi*(OD_1^4-ID_1^4)/32; % Cross-sectional polar moment of inertia


%% NiTi Tube 2 - Middle
OD_2 = 1.65e-3;
ID_2 = 1.40e-3;
len_2 = 80e-3;
len_cu_2 = 25e-3;
r2 = 85e-3;
u_y2_star = 1/r2;
I2 = pi*(OD_2^4-ID_2^4)/64;
J2 = pi*(OD_2^4-ID_2^4)/32;


%% NiTi Tube 3 - Innermost
OD_3 = 1.26e-3;
ID_3 = 1.1e-3;
len_3 = 120e-3;
len_cu_3 = 25e-3;
r3 = 65e-3;
u_y3_star = 1/r3;
I3 = pi*(OD_3^4-ID_3^4)/64;
J3 = pi*(OD_3^4-ID_3^4)/32;


%% Parameter array
len = [len_1,len_2,len_3];                  % Total length of NiTi tubes
len_cu = [len_cu_1,len_cu_2,len_cu_3];      % Curved length of NiTi tubes
ux_star = [0,0,0];                          % Pre-curvatures in the x-direction
uy_star = [u_y1_star,u_y2_star,u_y3_star];  % Pre-curvatures in the y-direction
I = [I1,I2,I3];                             % Cross-sectional moment of inertia
J = [J1,J2,J3];                             % Cross-sectional polar moment of inertia
