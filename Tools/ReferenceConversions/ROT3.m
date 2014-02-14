%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ROT3
% Date: 10/1/2013
%
% input:        theta = rotation angle [radians]
%
% output:       Rz = 3x3 rotation matrix about the Z (3rd) axis
%                   
% description:
% this fucntion computes the rotation matrix about the 3rd axis with
% respect to the passed rotation angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Rz] = ROT3(theta)

% Rotation about the third axis
Rz = [ cos(theta) sin(theta)  0;
      -sin(theta) cos(theta)  0;
      0           0           1];
 