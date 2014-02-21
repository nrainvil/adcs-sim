function [ b_v ] = rot_azel( x_v, z_v, az, el )
%AZEL_ROT create vector at Azimuth and Elevation
% az = aximuth (degrees) el = elevation (degrees)
    b_t_v = rot_srt(x_v, az, z_v);
    b_v = rot_srt(b_t_v, el, cross(b_t_v, z_v));
end

