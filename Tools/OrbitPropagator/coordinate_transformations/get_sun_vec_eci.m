function [sun_vec_eci] = get_sun_vec_eci(jd)

%% Constants
earth_eccentricity = 0.016734014;
earth_mean_motion_rdps = 1.99106168890789e-007;
earth_per_sun_vec = [ 0.226275654777816 -0.893696464855884 -0.387435100060598 ];
earth_per_sun_vec = earth_per_sun_vec/norm(earth_per_sun_vec);
earth_rot_axis = [ 0.0 -0.397756481161782 0.917491025402318 ];
earth_rot_axis = earth_rot_axis/norm(earth_rot_axis);
earth_perihelion_utc = 1294083118.8; % time of earth perihelion (UTC)
ecc_tolerance = 1.0e-10;
max_newton_iteration = 10;

%% Convert Julian date to unix time
[year month day hour min sec] = jd2date(jd);
t_utc = date2unixsecs(year,month,day,hour,min,sec);

%% Compute time since Earth was at perihelion
t_peri = t_utc - earth_perihelion_utc;

%% Compute mean anamoly
mean_anamoly = earth_mean_motion_rdps*t_peri;

%% Compute eccentric anamoly
% M = E - e*sin(E)
ecc0 = earth_eccentricity;
Etemp = mean_anamoly;
ratio = 1;
tol = ecc_tolerance;
n_iter = 0;
while abs(ratio) > tol && n_iter <= max_newton_iteration
    n_iter = n_iter+1;
    f_E = Etemp - ecc0*sin(Etemp) - mean_anamoly;
    f_Eprime = 1 - ecc0*cos(Etemp);
    ratio = f_E/f_Eprime;
    if abs(ratio) > tol
        Etemp = Etemp - ratio;
    else
        E = Etemp;
    end
end

true_anomaly = 2*atan2(sqrt(1+ecc0)*sin(E/2),sqrt(1-ecc0)*cos(E/2));

%% Rotate perihelion sun vector by true anomaly about the earth rotation axis
q = [earth_rot_axis*sin(true_anomaly/2) cos(true_anomaly/2)];
sun_vec_eci = Math_Q_Vec_Qinv_Mult(earth_per_sun_vec,q);
sun_vec_eci = sun_vec_eci/norm(sun_vec_eci);
