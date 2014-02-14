%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacob Cook
% 10/8/2013
% q2dcm.m
%
% input:    q = 4x1 quaternion
%
% output:   BI = 3x3 directin cosine matrix mapping ECI to body coordinates
%
% this function builds a DCM from the quaternion that maps a vector given
% in the ECI frame to the body frame.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BI] = q2dcm(q)

q1 = q(1,1);
q2 = q(2,1);
q3 = q(3,1);
q4 = q(4,1);

% Wertz p. 414
BI =[q1^2-q2^2-q3^2+q4^2      2*(q1*q2+q3*q4)      2*(q1*q3-q2*q4);
     2*(q1*q2-q3*q4)         -q1^2+q2^2-q3^2+q4^2  2*(q2*q3+q1*q4);
     2*(q1*q3+q2*q4)          2*(q2*q3-q1*q4)    -q1^2-q2^2+q3^2+q4^2];

