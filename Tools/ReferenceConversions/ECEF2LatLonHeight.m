%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ECEF2LatLonHeight
% Date: 10/1/2013
%
% input:        ECEF = 3x1 vector in ECEF frame [km]
%               R = Earth's radius [km]
%
% output:       phi = latitude [radians]
%               lambda = longitude [radians]
%               height = height from surface [km]
%                    
% description:
% This function computes the geocentric latitude longitude and height from
% the passed ECEF position coordinates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi, lambda, height] = ECEF2LatLonHeight(ECEF, R)

% Height
r = norm(ECEF); % magnitude of the position vector [km]
height = r-R; % height [km]
% Latitude
phi = asin(ECEF(3)/r); % latitude [radians]
% Longitude 
lambda_s = asin(ECEF(2)/(r*cos(phi))); 
lambda_c = acos(ECEF(1)/(r*cos(phi)));
lambda = quadrantCheck(lambda_s, lambda_c); % longitude [radians]
