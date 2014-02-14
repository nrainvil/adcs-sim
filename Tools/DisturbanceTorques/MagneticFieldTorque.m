%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 6/30/2013
% MagneticFieldTorque.m
%
% input:    cubeSat     = cubeSat struct containing geometric info
%           B_ECI       = 3x1 Magnetic Field vector expressed in the ECI 
%                         frame
%           q           = 4x1 quaternion defining the rotation between 
%                         from the body to the inertial frame.
%
% output:   N_mag       = 3x1 Magnetic Field torque in the body frame.
%
% This function computes the magnetic Field torque acting on the space 
% craft from its position in the ECI frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N_mag = MagneticFieldTorque(cubeSat, B_ECI, q)

% Express the B field in the body frame
[IB, BI] = q2dcm(q);
B_body = BI*B_ECI;
% compute the magnetic torque
N_mag = cross(cubeSat.MagneticDipole, B_body);
