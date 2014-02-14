%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 6/30/2013
% AttitudeProp.m
%
% input:    posECI_km   = satellite ECI position vector [km]
%           posECI_kmps = satellite ECI velocity vector [km/sec]
%
% output:   OI = 3x3 discrete cosine matrix to map ECI to orbital frame.
%           
% This function creates a DCM matrix that will map a ECI vector to Orbital
% coordinates where the axis are defined as follows:
%   X - along the velocity vector
%   Y - normal to the orbit
%   Z - completes the right hand rule (if circular it should point nadir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_IO] = ECI2Orbital(posECI_km, velECI_kmps)

% X-axis
i = velECI_kmps/norm(velECI_kmps); % along the velocity vector
% % Y-axis
orbitNormal = cross(posECI_km,velECI_kmps);
j = -orbitNormal/norm(orbitNormal); % orbit normal 
% Z-axis
k = cross(i,j);

% k = -posECI_km/norm(posECI_km);
% j = cross(k,i);
% DCM 
R_IO = [ i'; j'; k'];