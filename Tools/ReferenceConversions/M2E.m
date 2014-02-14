function [E] = M2E(M,e)
% p.73
tolerance = 1e-8;

if M > -pi && M < 0 || M > pi
    En = M-e;
else
    En = M+e;
end

Enp1 = 0; 
while abs(Enp1-En) > tolerance
    En = Enp1;
    Enp1 = En + (M-En+e*sin(En))/(1-e*cos(En));
end

E = Enp1;