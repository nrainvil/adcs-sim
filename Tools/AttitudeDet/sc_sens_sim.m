function [ B_mag_sens, I_sun_sens,  norm_sun_sens ] = sc_sens_sim(R_eci_body,  Sun_ECI, B_ECI, albedo_v)
%SC_SENSOR_SIM
%Simulate Attitude Determination Sensors
Sat_ECI_xhat = R_eci_body(:,1,:);
Sat_ECI_yhat = R_eci_body(:,2,:);
Sat_ECI_zhat = R_eci_body(:,3,:);

%3x Magnetometers
%B_mag_sens(1,:) = mag_sens(B_ECI,Sat_ECI_xhat);
%B_mag_sens(2,:) = mag_sens(B_ECI,Sat_ECI_yhat);
%B_mag_sens(3,:) = mag_sens(B_ECI,Sat_ECI_zhat);
B_mag_sens = mag_sens_3axis(B_ECI,Sat_ECI_xhat,Sat_ECI_zhat);

%14x Sun sensors (RAX2 positioning)
%norm_sun_sens_eci(1,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 17,-10);
%norm_sun_sens_eci(2,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,20);
%norm_sun_sens_eci(3,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -17,10);
%norm_sun_sens_eci(4,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -162,-10);
%norm_sun_sens_eci(5,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 180,20);
%norm_sun_sens_eci(6,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 162,-10);
%norm_sun_sens_eci(7,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 72,10);
%norm_sun_sens_eci(8,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 107,10);
%norm_sun_sens_eci(9,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 90,-20);
%norm_sun_sens_eci(10,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -107,10);
%norm_sun_sens_eci(11,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -72,10);
%norm_sun_sens_eci(12,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -90,-20);
%norm_sun_sens_eci(13,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,90);
%norm_sun_sens_eci(14,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,-90);
%norm_sun_sens_bf(1,:,:)  = rot_azel([1;0;0], [0;0;1], 17,-10);
%norm_sun_sens_bf(2,:,:)  = rot_azel([1;0;0], [0;0;1], 0,20);
%norm_sun_sens_bf(3,:,:)  = rot_azel([1;0;0], [0;0;1], -17,10);
%norm_sun_sens_bf(4,:,:)  = rot_azel([1;0;0], [0;0;1], -162,-10);
%norm_sun_sens_bf(5,:,:)  = rot_azel([1;0;0], [0;0;1], 180,20);
%norm_sun_sens_bf(6,:,:)  = rot_azel([1;0;0], [0;0;1], 162,-10);
%norm_sun_sens_bf(7,:,:)  = rot_azel([1;0;0], [0;0;1], 72,10);
%norm_sun_sens_bf(8,:,:)  = rot_azel([1;0;0], [0;0;1], 107,10);
%norm_sun_sens_bf(9,:,:)  = rot_azel([1;0;0], [0;0;1], 90,-20);
%norm_sun_sens_bf(10,:,:) = rot_azel([1;0;0], [0;0;1], -107,10);
%norm_sun_sens_bf(11,:,:) = rot_azel([1;0;0], [0;0;1], -72,10);
%norm_sun_sens_bf(12,:,:) = rot_azel([1;0;0], [0;0;1], -90,-20);
%norm_sun_sens_bf(13,:,:) = rot_azel([1;0;0], [0;0;1], 0,90);
%norm_sun_sens_bf(14,:,:) = rot_azel([1;0;0], [0;0;1], 0,-90);

%30 Deg pos
%+X
norm_sun_sens_eci(1,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 17,-10);
norm_sun_sens_eci(2,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,20);
norm_sun_sens_eci(3,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -17,10);
norm_sun_sens_bf(1,:,:)  = rot_azel([1;0;0], [0;0;1], 17,-10);
norm_sun_sens_bf(2,:,:)  = rot_azel([1;0;0], [0;0;1], 0,20);
norm_sun_sens_bf(3,:,:)  = rot_azel([1;0;0], [0;0;1], -17,10);
%-X
norm_sun_sens_eci(4,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -162,-10);
norm_sun_sens_eci(5,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 180,20);
norm_sun_sens_eci(6,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 162,-10);
norm_sun_sens_bf(4,:,:)  = rot_azel([1;0;0], [0;0;1], -162,-10);
norm_sun_sens_bf(5,:,:)  = rot_azel([1;0;0], [0;0;1], 180,20);
norm_sun_sens_bf(6,:,:)  = rot_azel([1;0;0], [0;0;1], 162,-10);
%-Y
norm_sun_sens_eci(7,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 90,-30);
norm_sun_sens_eci(8,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 90,30);
norm_sun_sens_eci(9,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 120,0);
norm_sun_sens_bf(7,:,:)  = rot_azel([1;0;0], [0;0;1], 90,-30);
norm_sun_sens_bf(8,:,:)  = rot_azel([1;0;0], [0;0;1], 90,30);
norm_sun_sens_bf(9,:,:)  = rot_azel([1;0;0], [0;0;1], 120,0);
%+Y
norm_sun_sens_eci(10,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -90,-30);
norm_sun_sens_eci(11,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -90,30);
norm_sun_sens_eci(12,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -120,0);
norm_sun_sens_bf(10,:,:) = rot_azel([1;0;0], [0;0;1], -90,-30);
norm_sun_sens_bf(11,:,:) = rot_azel([1;0;0], [0;0;1], -90,30);
norm_sun_sens_bf(12,:,:) = rot_azel([1;0;0], [0;0;1], -120,0);
%+Z
norm_sun_sens_eci(13,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,90);
norm_sun_sens_bf(13,:,:) = rot_azel([1;0;0], [0;0;1], 0,90);
%-Z
%N/A

sim_length = length(Sat_ECI_xhat(1,:));


[I_sun_sens, I_sun_sens_alpha]   = sun_sens(Sun_ECI,albedo_v,norm_sun_sens_eci); 
norm_sun_sens = norm_sun_sens_bf;
end
