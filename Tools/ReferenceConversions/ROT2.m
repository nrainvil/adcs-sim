%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ROT2
% Date: 10/1/2013
%
% input:        theta = rotation angle [radians]
%
% output:       Ry = 3x3 rotation matrix about the Y (2nd) axis
%                   
% description:
% this fucntion computes the rotation matrix about the 2nd axis with
% respect to the passed rotation angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ry] = ROT2(theta)

% Computes DCM for rotation about the Y axis
Ry = [cos(theta) 0 -sin(theta); 
      0          1  0; 
      sin(theta) 0  cos(theta)];