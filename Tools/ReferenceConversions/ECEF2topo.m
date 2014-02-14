%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ECEF2topo
% Date: 10/1/2013
%
% input:        ECEF = 3x1 position of satellite in the ECEF frame
%               phi = ground site latitude [radians]
%               lambda = ground site longitude [radians]
%               height = ground site height from surface [km]
%               R = Earth's radius [km]
%
% output:       r_topo = 3x1 vector containing range, elevation and azimuth
%                    
% description:
% This function computes the topographical range, elevation and azimuth of
% a satellite relative to a ground site, gine the position of the satellite
% in ECEF coordinates and the latitude longitude and height of the ground
% station.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_topo] = ECEF2topo(r_ECEF, phi, lambda, height, R)

% Get ECEF coordinates of the gorund site
r_site_ECEF = LatLonHeight2ECEF(phi, lambda, height, R); % [km]

% compute slant-range vector from site to satellite
rho_ECEF = r_ECEF - r_site_ECEF; % [km]

% rotate to SEZ frame 
rho_SEZ = ROT2(pi-phi)*ROT3(lambda)*rho_ECEF; % [km]

% Range
rho = norm(rho_SEZ); % [km]

% Elevation
el_s = asin(rho_SEZ(3)/rho); 
el_c = acos((rho_SEZ(1)^2+rho_SEZ(2)^2)^.5/rho);
el = quadrantCheck(el_s, el_c); % [radians]

% Azimuth
beta_s = asin(rho_SEZ(2)/(rho_SEZ(1)^2+rho_SEZ(2)^2)^.5);
beta_c = acos(-rho_SEZ(1)/(rho_SEZ(1)^2+rho_SEZ(2)^2)^.5);
beta = quadrantCheck(beta_s, beta_c); % [radians]

% group and ship
r_topo = [rho el beta]';

