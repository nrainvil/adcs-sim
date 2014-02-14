%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ECEF2ECI
% Date: 10/1/2013
%
% input:        ECEF = 3x1 vector in ECEF frame [km]
%               theta_GST = Greenwich Sidereal Time [radians]
%
% output:       ECI = 3x1 vector in ECI frame [km]
%                    
% 
% description:
% this function computes the ECI vector given an ECEF vector and the
% Greenwich Sidereal Time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_ECI] = ECEF2ECI(r_ECEF, theta_GST)

% Rotate about the Z axis
r_ECI = ROT3(-theta_GST)*r_ECEF;