function [C] = ROT313(RAAN, omega, i)

cR = cos(RAAN);
sR = sin(RAAN);
cw = cos(omega);
sw = sin(omega);
ci = cos(i);
si = sin(i);

C = [cR*cw-sR*sw*ci -cR*sw-sR*cw*ci  sR*si;
     sR*cw+cR*sw*ci -sR*sw+cR*cw*ci -cR*si;
     sw*si           cw*si           ci];
 