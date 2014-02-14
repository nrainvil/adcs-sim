%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 6/30/2013
% AerodynamicTorque.m
%
% input:    cubeSat     = cubeSat struct containing geometric info
%           velECI_kmps = 3x1 velocity vector in the ECI frame
%           posECEF_km  = 3x1 position vector in the ECEF frame
%           Time        = datenum serial time 
%           q           = 4x1 quaternion defining the rotation between from
%                         the body to the inertial frame.
%
% output:   N_aero      = 3x1 aerodynamic drag torque in the body frame.
%
% This function computes the aerodynamic drag torque acting on the space 
% craft from its position, velocity and geometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N_aero = AerodynamicTorque(cubeSat, velECI_kmps, posECEF_km, time, q)

% Convert Inertial Velocity vector to the body frame
[IB, BI] = q2dcm(q);
V_b = BI*velECI_kmps;
V_mag = norm(V_b);
V_hat = V_b/V_mag;
V_mag = V_mag*1000; % convert to meters/sec

% Time 
[year, month, date, hour, min, sec] = datevec(time);
UTsec = hour*3600 + min*60 + sec;
dayOfYear = time - datenum(year, 1, 1, hour, min, sec)+1;

% Atmospheric Density
posECEF_m = posECEF_km*1000; % convert to meters
LatLonAlt = ecef2lla(posECEF_m'); % get Lat Lon and Altitude from ECEF frame
Lat_m = LatLonAlt(1)*111.32e3; % m
Lon_deg = LatLonAlt(2); % degrees
Alt_m = LatLonAlt(3); % m
% use matlab nrlmsise00 to get atmoshperic density
[T, rhoVec] = atmosnrlmsise00(Alt_m, Lat_m, Lon_deg, year, dayOfYear, UTsec,'none');
rho = rhoVec(6);

%Drag Coefficient
C_d = cubeSat.DragCoefficient;

% Book keeping
nSides = length(cubeSat.SurfaceArea);
N_a = zeros(3,nSides);

% Compute drag torques produced by each side
for i = 1:nSides
    % get side specific info
    A = cubeSat.SurfaceArea(i); % surface area of the side
    N_hat = cubeSat.SurfaceNormal(:,i); % Normal vector in body coordinate
    R_cp = cubeSat.CenterOfPressure(:,i); % Center of pressure in body coordinates
    
    if N_hat'*V_hat <= 0 % face not flying into the wind
        N_a(:,i) = 0;
    else
        % Compute drag force
        F_a = -0.5*C_d*rho*V_mag^2*A*(N_hat'*V_hat)*V_hat;
        % Compute drag torque
        N_a(:,i) = cross(R_cp,F_a);
    end
end

% Sum drag torques
N_aero = sum(N_a,2);