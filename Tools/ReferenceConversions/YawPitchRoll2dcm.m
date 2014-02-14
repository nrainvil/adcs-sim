function [dcm] = YawPitchRoll2dcm(ypr)

cr = cos(ypr(3));
sr = sin(ypr(3));
cp = cos(ypr(2));
sp = sin(ypr(2));
cy = cos(ypr(1));
sy = sin(ypr(1));

dcm = [cp*cy           cp*sy          -sp;
       sr*sp*cy-cr*sy  sr*sp*sy+cr*cy  sr*cp;
       cr*sp*cy+sr*sy  cr*sp*sy-sr*cy  cr*cp];