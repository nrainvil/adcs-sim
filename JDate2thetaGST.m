function [ theta01 ] = JDate2thetaGST( Jdate )
%JDATE2THETAGST - John Stark
T = (Jdate - 2451545)/36525; % Find T_UT1

% Calculate theta_GMST from the equation given in Vallado
theta01 = 67310.54841+(876600+8640184.812866)*T+0.093104*T^2-6.2e-6*T^3;
while theta01 > 86400
    theta01 = theta01 - 86400; % Correct such that the value is in the 
                               % first rotation. Quantity is in seconds
end
theta01 = theta01/240*pi/180; % Convert to radians


end

