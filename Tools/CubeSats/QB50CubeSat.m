function [CubeSat] = QB50CubeSat()

% %% Geometric and Mechanical Parameters of CubeSat
CubeSat.I_body = [   0.0033   0    0; 
                        0   0.0101 0; 
                        0   0    0.0101]; % kg*m^2
CubeSat.mass = 2; % kg
CubeSat.Dimensions = [.225 .100 .100]'; % dimension [x y z] in meters
CubeSat.GeometricCenter = [0 0 0]'; % geometric center location in body coordinates
CubeSat.CG = [0.0 0.0 0.0 ]'; % center of gravity location in body coordinates

%%% Center of pressure %%%
% here the center of pressure is defined 
CubeSat.CenterOfPressure = [ 0.1125  0     0;       % X+
                            -0.1125  0     0;       % X-
                             0       0.05  0;       % Y+
                             0      -0.05  0;       % Y-
                             0       0     0.05;    % Z+
                             0       0    -0.05]';  % Z-

% Surface Area           X+    X-     Y+      Y-      Z+      Z-
CubeSat.SurfaceArea = [ 0.01, 0.01, 0.0225, 0.0225, 0.0225, 0.0225]; 
CubeSat.SurfaceNormal = [ 1  0  0;   % X+
                         -1  0  0;   % X-
                          0  1  0;   % Y+
                          0 -1  0;   % Y-
                          0  0  1;   % Z+
                          0  0 -1]'; % Z-

CubeSat.DragCoefficient = 2.5;
CubeSat.MagneticDipole = [0.001 0 0]'; % A*m^2, Magnetic dipole along x axis

%% Reaction Wheel Parameters
CubeSat.RWMass = 0.150; % kg, mass of reaction wheel
CubeSat.RWRadius = 0.0215; % m, radius of reaction wheel
CubeSat.RWHeight = 0.018; % m, height of reaction wheel
I_rw_x = 1/12*CubeSat.RWMass*(3*CubeSat.RWRadius^2+CubeSat.RWHeight^2);
I_rw_y = CubeSat.RWMass*CubeSat.RWRadius^2/2;
I_rw_z = I_rw_x;
CubeSat.RWInertiaMatrix = [I_rw_x 0      0;  % kg*m^2, reaction wheel Inertia Tensor
                           0      I_rw_y 0;  % in body coordinates
                           0      0      I_rw_z];
CubeSat.RWMomentumVector = [0 1 0]'; % in body coordinates

%% Control Paramters
CubeSat.K_pitch = 0; %.0001; % pitch gain
CubeSat.rho_pitch = 0; % damping ratio of control loop
CubeSat.PitchTimeConstant = 2*CubeSat.rho_pitch*(CubeSat.K_pitch/I_rw_y)^.5;
                           
