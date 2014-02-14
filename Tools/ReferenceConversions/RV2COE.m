%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 6/30/2013
% RV2COE.m
%
% input:    posECI_km   = satellite ECI position vector [km]
%           posECI_kmps = satellite ECI velocity vector [km/sec]
%
% output:   COE         = Stuct containing the classical orbital elements
%                           a     - semi-major axis
%                           e     - eccentricity
%                           i     - inclination
%                           omega - argument of perigee
%                           RAAN  - right ascension of the ascending node
%                           nu    - true anomaly
%           
% This function computes the classical orbital elements from the ECI
% position and velocity given in km an kmps respectively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [COE] = RV2COE(posECI_km,velECI_kmps)

global MU_EARTH

% vector magnitudes
r_mag = (posECI_km'*posECI_km)^.5; % km
v_mag = (velECI_kmps'*velECI_kmps)^.5; % km/s
v_r = (velECI_kmps'*posECI_km/r_mag); % km/s

% angular momentum
h_vec = cross(posECI_km,velECI_kmps);
h = (h_vec'*h_vec)^.5;

% semi major axis
xi = v_mag^2/2 - MU_EARTH/r_mag;
a = -MU_EARTH/(2*xi);

% inclination
i = acos(h_vec(3,1)/h);

% Node line
N = cross([0 0 1]', h_vec);
N_mag = (N'*N)^.5;

% right ascension of the ascending node
RAAN = acos(N(1,1)/N_mag);
if N(2,1) < 0
    RAAN = 2*pi - RAAN;
end

% eccentricity
e_vac = 1/MU_EARTH*((v_mag^2-MU_EARTH/r_mag)*posECI_km - r_mag*v_r*velECI_kmps);
e = (e_vac'*e_vac)^.5;

% argument of perigee
omega = acos((N/N_mag)'*(e_vac/e));
if e_vac(3,1) < 0 
    omega = 2*pi-omega;
end

% true anomoly
nu = acos((e_vac/e)'*(posECI_km/r_mag));
if v_r< 0
    nu = 2*pi - nu;
end

% Pack in a structure 
COE.a = a;
COE.e = e;
COE.i = i;
COE.omega = omega;
COE.RAAN = RAAN;
COE.nu = nu;

