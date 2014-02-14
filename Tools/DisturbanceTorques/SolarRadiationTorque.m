%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 6/30/2013
% SolarRadiationTorque.m
%
% input:    cubeSat     = cubeSat struct containing geometric info
%           Sun_ECI     = 3x1 sun vector expressed in the ECI frame
%           q           = 4x1 quaternion defining the rotation between 
%                         from the body to the inertial frame.
%
% output:   N_solar     = 3x1 solar radiation torque in the body frame.
%
% This function computes the solar radiation torque acting on the space 
% craft from its orientation relative to the sun and geometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N_solar = SolarRadiationTorque(cubeSat, Sun_ECI, q)

% Sun vector in body coordinates
[IB, BI] = q2dcm(q);
S_hat = BI*Sun_ECI;

% Solar Pressure constant
P = 4.4e-6; % kg/(m*s) mean momentum flux
C_s = 0.4;
C_d = 0.2;

% Book keeping
nSides = length(cubeSat.SurfaceArea);
N_s = zeros(3,nSides);

% Compute torques produced by each side
for i = 1:nSides
    N_hat = cubeSat.SurfaceNormal(:,i);
    A = cubeSat.SurfaceArea(i);
    R_cp = cubeSat.CenterOfPressure(:,i);
   % check for the 
   if norm(S_hat) == 0
       N_s(:,:) = 0;
       break
   else
       theta = acos(S_hat'*N_hat);
       if theta < 0 || theta > pi/2 % surface not illuminated
           F_i(1:3,1) = 0;
       else
           % Solar radiation force
           F_i = -P*A*cos(theta)*((1-C_s)*S_hat + 2*(C_s*cos(theta)+1/3*C_d)*N_hat);
       end
       % solar radiation torque
       N_s(:,i) = cross(R_cp, F_i);
   end
end

% sum of all torques
N_solar = sum(N_s,2);




