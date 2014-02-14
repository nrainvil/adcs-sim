function [R V] = COE2RV(OE)

mu = 3.986005e5; % km^3/s^2 Geocentric Gravitational Constant
p = OE.a*(1 - OE.e^2);

Rpqw = [p*cos(OE.nu)/(1+OE.e*cos(OE.nu)) p*sin(OE.nu)/(1+OE.e*cos(OE.nu)) 0]';

Vpqw = [-(mu/p)^.5*sin(OE.nu) (mu/p)^.5*(OE.e+cos(OE.nu)) 0]';

IP = ROT313(OE.RAAN, OE.omega, OE.i);

R = IP*Rpqw;
V = IP*Vpqw;