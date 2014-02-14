function [x_km_mag] = ecf2mag(x_km_ecf,mag_parm)

% earth_dipole_ecf = [0.06859, -0.18602, 0.98015];
% earth_dipole_ecf = earth_dipole_ecf/norm(earth_dipole_ecf);
earth_dipole_ecf = mag_parm.dipole_axis_unit_ecf;
dipole_offset_km_ecf = mag_parm.dipole_pos_Re_ecf*mag_parm.Re_km;

%% Subtract dipole position offset
x_tmp = x_km_ecf - dipole_offset_km_ecf;

%% Find the rotation axis between ECF Z and the dipole axis
v = cross([0 0 1],earth_dipole_ecf);
v = v/norm(v); % this is the MAG Y axis

%% Find the dipole tilt angle
theta_tilt = acos(dot([0 0 1],earth_dipole_ecf));

%% Rotate ECF Y to MAG Y about Z
theta_z = acos(dot([0 1 0],v));
evec = cross([0 1 0],v);
evec = evec/norm(evec);

qz = [evec*sin(theta_z/2) cos(theta_z/2)];

%% Rotate about MAG Y to align Z with the dipole axis
qy = [[0 1 0]*sin(theta_tilt/2) cos(theta_tilt/2)];

q = Math_Q_Q_Mult(qz,qy);

x_km_mag = Math_Qinv_Vec_Q_Mult(x_tmp,q);