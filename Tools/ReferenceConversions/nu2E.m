function [E] = nu2E(nu,e)
% p. 85
Es  = asin(sin(nu)*(1-e^2)^.5/(1+e*cos(nu)));
Ec = acos((e+cos(nu))/(1+e*cos(nu)));

if Es < 0 && Ec < pi/2      % quadrant IV
    E = 2*pi - Ec;
elseif Es < 0 && Ec > pi/2  % quadrant III
    E = pi - Es;
elseif Es > 0 && Ec > pi/2  % quadrant II
    E = Ec;
else                        % quadrant I
    E = Es;
end
