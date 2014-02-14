%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ECI2ECEF
% Date: 10/1/2013
%
% input:        ECI = 3x1 vector in ECI frame [km]
%               theta_GST = Greenwich Sidereal Time [radians]
%
% output:       ECEF = 3x1 vector in ECEF frame [km]
%                    
% 
% description:
% this function computes the ECEF vector given an ECI vector and the
% Greenwich Sidereal Time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ECEF] = ECI2ECEF(ECI, theta_GST)

% Rotate about the Z-axis 
ECEF = ROT3(theta_GST)*ECI;