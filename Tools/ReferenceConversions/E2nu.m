function [nu] = E2nu(E, e)
% p. 85
nus = asin(sin(E)*(1-e^2)^.5/(1-e*cos(E)));
nuc = acos((cos(E)-e)/(1-e*cos(E)));

if nus < 0 && nuc < pi/2 % quadrant IV
    nu = 2*pi - nuc;
elseif nus < 0 && nuc > pi/2 % quadrant III
    nu = pi - nus;
elseif nus > 0 && nuc > pi/2 % quadrant II
    nu = nuc;
else
    nu = nus;
end
