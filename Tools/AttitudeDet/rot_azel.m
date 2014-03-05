function [ b_v ] = rot_azel( x_v, z_v, az, el )
% AZEL_ROT create vector at Azimuth and Elevation
% az = aximuth (degrees) el = elevation (degrees)
    for i = 1:length(x_v(1,:))
        b_t_v(:,i) = rot_srt(x_v(:,i), z_v(:,i), az);
        b_v(:,i) = rot_srt(b_t_v(:,i), cross(b_t_v(:,i), z_v(:,i)), el);
    end
end

