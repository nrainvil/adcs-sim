%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: quadrantCheck
% Date: 10/1/2013
%
% input:        sa = arcsine value [radians]
%               ca = arccosine value [radians]
%
% output:       a = final angle [radians]
%                   
% description:
% This function computes the true 4 quadrant angle from the passed inverse
% sine and cosine values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a] = quadrantCheck(sa, ca)

% Check Quadrant
if sa < 0 && ca < pi/2      % quadrant IV
    a = 2*pi - ca;
elseif sa < 0 && ca > pi/2  % quadrant III
    a = pi - sa;
elseif sa > 0 && ca > pi/2  % quadrant II
    a = ca;
else                        % quadrant I
    a = sa;
end
