%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 10/8/2013
% OrbitAngularVelocity.m
%
% input:    posECI_km   = satellite ECI position vector [km]
%           posECI_kmps = satellite ECI velocity vector [km/sec]
%
% output:   OMEGA = instantaneous angular rate of satellites position
%
% This function computes the instananeous angular velocity of the satallite
% given the current position and velocity of the spacecraft.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [OMEGA] = OrbitAngularVelocity(posECI_km, velECI_kmps)

global MU_EARTH

% Get classical orbital elements from state vector
OE = RV2COE(posECI_km,velECI_kmps);
r_mag = (posECI_km'*posECI_km)^.5;

% semi-minor axis
b = OE.a*(1-OE.e^2 )^.5;
% orbital period
T = 2*pi/MU_EARTH^.5*OE.a^(3/2);
% mean motion
n = 2*pi/T;
% Instantaneous Angular rate of orbit tangent
OMEGA = n*OE.a*b/r_mag^2; % [rad/sec]