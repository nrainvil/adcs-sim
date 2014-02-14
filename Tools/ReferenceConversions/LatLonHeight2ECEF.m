%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: LatLonHeight2ECEF
% Date: 10/1/2013
%
% input:        phi = latitude [radians]
%               lambda = longitude [radians]
%               height = height from surface [km]
%               R = Earth's radius [km]
%
% output:       ECEF = 3x1 vector in ECEF frame [km]
%                    
% description:
% This function computes the ECEF position coordinates from the geocentric 
% latitude longitude and height.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_ECEF] = LatLonHeight2ECEF(phi, lambda, height, R)

% Compute position radius
r = R + height; % magnitude of position vector [km]
% Compute ECEF Coordinates
r_ECEF(1,1) = r*cos(phi)*cos(lambda); % X ECEF [km]
r_ECEF(2,1) = r*cos(phi)*sin(lambda); % Y ECEF [km]
r_ECEF(3,1) = r*sin(phi); % Z ECEF [km]