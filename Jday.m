function [JD] = Jday(Y,M,D,UT)
%Convert Y/M/D to J2000

JD0 = 367.*Y-floor(7/4*(Y+floor((M+9)/12)))+floor(275*M/9)+D+1721013.5;
JD = JD0+UT/24;

end

