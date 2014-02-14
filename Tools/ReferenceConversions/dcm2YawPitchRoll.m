function [YawPitchRoll] = dcm2YawPitchRoll(dcm)

yaw = atan2(dcm(1,2),dcm(1,1));
pitch = asin(-dcm(1,3));
roll = atan2(dcm(2,3),dcm(3,3));
YawPitchRoll = [yaw pitch roll]';