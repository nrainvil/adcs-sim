function Ng = GravityGradientTorque(I_body, posECEF_km, q)


%% Global Variable
global MU_EARTH

[IB, BI] = q2dcm(q);
pos_b = BI*posECEF_km;
R_mag = norm(pos_b); % magnitude of position vector
R_hat = pos_b/R_mag; % position unit vector

R_mag_m = R_mag*1000; % put radius vector in meters

Ng = 3*MU_EARTH/R_mag_m^3*cross(R_hat, I_body*R_hat);


