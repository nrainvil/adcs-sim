%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 10/5/2013
% OE2TLE.m
%
% input:    time = [1xN double] (UTC) time in matlab datenum format
%           TLE  =  [1x1 structure] contains TLE information. Has fields:
%                   header: [char] indentifying info for the TLE. Not used in this function.
%                   line1: [1x69 char] line 1 of the TLE
%                   line2: [1x69 char] line 2 of the TLE
%
% output:   Orbit = struct containing relevant orbital information
%
% This function is a simple wrapper for David Gerhardts Orbital propagator
% to throw the propagators outputs into a structure that is used within the
% attitude propagator.  This helps to easy the passing of variables between
% function. TODO implement a specific Orbit class that can specify return
% much more from the info given below such as the true anomoly orbit frame
% conversions and the like.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Orbit ] = GenOrbit(time, TLE)

% use David Gerhardts Orbit propagator to generate Orbit 
[posECI_km, posECEF_km, velECI_kmps, S_ECI, B_ECI]=GenPosVelSunMag(time,TLE);

% put orbit parameters in Orbit Structure
Orbit.posECEF_km = posECEF_km;
Orbit.posECI_km = posECI_km;
Orbit.velECI_kmps = velECI_kmps;
Orbit.Sun_ECI = S_ECI;
Orbit.B_ECI = B_ECI;
Orbit.Time = (time - time(1))*24*3600; % [sec]

end

