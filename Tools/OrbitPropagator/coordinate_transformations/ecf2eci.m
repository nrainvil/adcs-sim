function x_eci = ecf2eci(x_ecf,jd)

%% Inputs
% x_ecf: vector in ECF frame (not velocity or acceleration, rotating frame)
% jd: Julian date

%% Outputs
% x_eci: vector in ECI frame

%% Notes
% See ECI2ECEF, results are slightly different, but agree to within 1e-5
% deg
% Extra steps required if converting a velocity or acceleration because ECF
% is a rotating frame

%% Constants
ecf2eci_epoch_time = 1316563312.8565;   % unix time (UTC)
% ecf2eci_epoch_angle = 0.0;              % rad
earth_rot_rate = 7.2921158220e-005;     % rad/sec

%% Convert Julian date to UTC
[year month day hour min sec] = jd2date(jd);
t_utc = date2unixsecs(year,month,day,hour,min,sec);

%% Compute angle between ECI and ECF
theta_eci2ecf_rad = (t_utc - ecf2eci_epoch_time)*earth_rot_rate;

q_eci2ecf = [[0 0 1]*sin(theta_eci2ecf_rad/2) cos(theta_eci2ecf_rad/2)];

%% Rotate vector about the Z axis
x_eci = Math_Q_Vec_Qinv_Mult(x_ecf,q_eci2ecf);