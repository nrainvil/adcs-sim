function [CubeSat] = BaseCubeSat()

%% %%%%%%%%%%%%%%%%%%%%%%Dimensions and Mass %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CubeSat.mass = 2;                       % kg
CubeSat.Dimensions = [.225 .100 .100]'; % dimension [x y z] in meters
CubeSat.GeometricCenter = [0 0 0]';     % geometric center location in body coordinates
CubeSat.CG = [0.0 0.0 0.0 ]';           % center of gravity location in body coordinates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%% Inertia Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is estimated by assuming the cubeSat is a solid block and the mass
% is evenly distributed about the volume

CubeSat.I_body = [   0.0033   0    0; 
                        0   0.0101 0; 
                        0   0    0.0101]; % kg*m^2
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
%% %%%%%%%%%%%%%%%%%%%% Center of pressure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the center of pressure is defined for each side of the cubesat.

CubeSat.CenterOfPressure = [ 0.1125  0     0;       % X+
                            -0.1125  0     0;       % X-
                             0       0.05  0;       % Y+
                             0      -0.05  0;       % Y-
                             0       0     0.05;    % Z+
                             0       0    -0.05]';  % Z-
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%% Surface Area and Surface Normal %%%%%%%%%%%%%%%%%%%%%%  
% The surface area parameter gives the surface area of each side of the
% cubesat and the Surface normal 

% Surface Area           X+    X-     Y+      Y-      Z+      Z-
CubeSat.SurfaceArea = [ 0.01, 0.01, 0.0225, 0.0225, 0.0225, 0.0225]; 

CubeSat.SurfaceNormal = [ 1  0  0;   % X+
                         -1  0  0;   % X-
                          0  1  0;   % Y+
                          0 -1  0;   % Y-
                          0  0  1;   % Z+
                          0  0 -1]'; % Z 
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Aditional Paramaters %%%%%%%%%%%%%%%%%%%%%%%%
CubeSat.DragCoefficient = 2.5;
CubeSat.MagneticDipole = [0.001 0 0]'; % A*m^2, Magnetic dipole along x axis

                           
