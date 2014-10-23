%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nicholas Rainville, Jacob Cook
% 10/9/2013
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all;
%close all;
warning('off','Control:analysis:LsimStartTime');

%% Load Orbit and attitude toolbox
wdir = pwd;
toolsDir = strcat(wdir,'/Tools');
toolPath = genpath(toolsDir);
addpath(toolPath)
cd(wdir)

clear baseDir toolPath toolsDir wdir

%% Global Variables
global R_EARTH MU_EARTH
R_EARTH = 6378.1363;        % [km]
MU_EARTH = 398600.4415;     % [km^3/s^2]
AU = 149597871;             % [km]

%% Set Classical Orbital Elements
% Note: We start at closest point to south pole to line up with our K
OE.a        = 370 + R_EARTH;    % semi-major axis [km]
OE.e        = 0;                % eccentricity
OE.i        = deg2rad(90);      % inclination [rad]
OE.omega    = deg2rad(0);       % argument of perigee [rad]
OE.RAAN     = deg2rad(90); %deg2rad(30);       % right ascension of the ascending node [rad]
OE.nu       = deg2rad(-90);     % true anomoly [rad]

%% TLE
SatNum = '02059';       % Satellite Number
SatClass = 'U';         % Satellite Classification
IntDsgntr = '12001A  '; % International Designator
Epoch = datenum('01 May 2015 00:00:00'); % [UTC] 
TLE = BuildTLE(SatNum,SatClass,IntDsgntr,Epoch,OE);

%% Simulation time
julienDate = Jday(2015,5,1,0);
j2000_offset = 2451545;
j2000Date = julienDate -j2000_offset;
startTime = datenum('01 May 2015 00:00:00');
stopTime  = datenum('01 May 2015 01:00:00');
%stopTime  = datenum('01 May 2015 01:32:00');
timeStep = .5;%.2; % [sec] (5 Hz)
time = startTime:timeStep/(3600*24):stopTime;

%Convert to GST using constants from http://aa.usno.navy.mil/faq/docs/GAST.php
time_hours = (time-floor(startTime))*24;
time_GST = 6.697374558 + 0.06570982441908*j2000Date + 1.00273790935*time_hours;
theta_GST = ((time_GST)/24)*(2*pi()); 

Orbit = GenOrbit(time, TLE);
clear SatNum SatClass IntDsgntr Epoch time

%% ECI to Satellite Body Frame DCM
Sat_ECI_xhat = zeros(size(Orbit.posECI_km));
Sat_ECI_zhat = zeros(size(Orbit.posECI_km));
Sat_ECI_yhat = zeros(size(Orbit.posECI_km));

Sat_ECI_zhat(1,:) = -1*Orbit.velECI_kmps(1,:)./sqrt(sum(abs(Orbit.velECI_kmps).^2,1));
Sat_ECI_zhat(2,:) = -1*Orbit.velECI_kmps(2,:)./sqrt(sum(abs(Orbit.velECI_kmps).^2,1));
Sat_ECI_zhat(3,:) = -1*Orbit.velECI_kmps(3,:)./sqrt(sum(abs(Orbit.velECI_kmps).^2,1));

for i=1:length(Orbit.velECI_kmps(1,:))
    Sat_ECI_xhat(:,i) = ortho_norm(-Orbit.posECI_km(:,i),Sat_ECI_zhat(:,i)); % Nadir pointing
    %Sat_ECI_xhat(:,i) = ortho_norm(Orbit.Sun_ECI_noecl(:,i),Sat_ECI_zhat(:,i)); % Sun Pointing
    Sat_ECI_yhat(:,i) = cross(Sat_ECI_zhat(:,i),Sat_ECI_xhat(:,i));
end

R_eci_body(:,1,:) = Sat_ECI_xhat;
R_eci_body(:,2,:) = Sat_ECI_yhat;
R_eci_body(:,3,:) = Sat_ECI_zhat;

for i=1:length(R_eci_body(1,1,:))
	S_true_bf(:,i) = reshape(R_eci_body(:,:,i),3,3)'*Orbit.Sun_ECI(:,i);
	B_true_bf(:,i) = reshape(R_eci_body(:,:,i),3,3)'*Orbit.B_ECI(:,i);
end

%% ECI Body Rates - Radians/sec
for i=2:length(Sat_ECI_xhat)
    R_m = reshape(R_eci_body(:,:,i-1),3,3);
    R_p = reshape(R_eci_body(:,:,i),3,3);
    dr_dt = ((R_p-R_m)/timeStep); %Body Frame Rates
    w_i = dr_dt*R_p'; %ECI Frame
    Sat_ECI_xrate(:,i) = w_i(3,2);
    Sat_ECI_yrate(:,i) = w_i(1,3);
    Sat_ECI_zrate(:,i) = w_i(2,1);
    R_dot_eci_body(:,:,i) = w_i;
end
Sat_ECI_xrate(:,1) = Sat_ECI_xrate(:,2);
Sat_ECI_yrate(:,1) = Sat_ECI_yrate(:,2);
Sat_ECI_zrate(:,1) = Sat_ECI_zrate(:,2);

R_dot_eci_body(:,:,1) = reshape(R_dot_eci_body(:,:,2),3,3);

%% Sensor Simulation
%Continuous
sim_length = length(Orbit.Sun_ECI);
load('mean_val.mat');
fprintf(1,'Sensor Sim:      ');
for i=1:sim_length
	Sun_ECEF = ECI2ECEF(Orbit.Sun_ECI(:,i),theta_GST(i));
	fprintf(1,'\b\b A');
        albedo_map = albedo(Orbit.posECEF_km(:,i).*1000,Sun_ECEF*AU*1000,mean_val);
	albedo_v = albedo_vec(albedo_map,Orbit.posECEF_km(:,i).*1000);
	fprintf(1,'\b\b S');
        [B_mag_sens_single, I_sun_sens_single,  norm_sun_sens_single] = sc_sens_sim(R_eci_body(:,:,i), Orbit.Sun_ECI(:,i),Orbit.B_ECI(:,i),albedo_v);
	B_mag_sens(:,i) = B_mag_sens_single;
        I_sun_sens(:,i) = I_sun_sens_single;
        norm_sun_sens = norm_sun_sens_single;
	fprintf(1,'\b\b\b\b\b%02d%% E',floor((i/sim_length)*100));
end
fprintf(1,'\n');

%Rate gyros - lsim requires time history
[G_rate_sens] = sc_rate_sim(R_dot_eci_body);

%Discrete
fs_sens = 2; %Hz
[B_mag_sens_rs, I_sun_sens_rs, G_rate_sens_rs, B_ECI_rs, Sun_ECI_rs, Time_rs] = adc_sim(fs_sens, B_mag_sens, I_sun_sens, G_rate_sens, Orbit);
sens_length = length(Time_rs);

%% Attitude Estimate
fprintf(1,'Attitude Estimate:    \n');
%Estimate sun vector
I_sun_sens_ds = reshape(I_sun_sens_rs.data,length(I_sun_sens(:,1)),sens_length);
%[S_sens_est_bf, S_sens_num_bf] = est_sun_sens(I_sun_sens_ds, norm_sun_sens);

%Estimate magnetic field vector
B_mag_sens_ds = reshape(B_mag_sens_rs.data,length(B_mag_sens(:,1)),sens_length);
B_sens_est = B_mag_sens_ds;

%Estimate body rates from gyros
G_rate_sens_ds = reshape(G_rate_sens_rs.data,length(G_rate_sens(:,1)),sens_length);
G_rate_est = G_rate_sens_ds;

% Generate Cosine Matrix
B_ECI_ds = reshape(B_ECI_rs.data,length(Orbit.B_ECI(:,1)),sens_length);
Sun_ECI_ds = reshape(Sun_ECI_rs.data,length(Orbit.Sun_ECI(:,1)),sens_length);
R_test(:,:,1) = R_eci_body(:,:,1); %START FROM KNOWN STATE (REMOVE)
timestep = 1/fs_sens;
drift_bias = [0;0;0];
z_k = 0;
P_k = eye(4,4);
for k=1:sens_length    
    %EKF
%    b_eci_k = B_ECI_ds(:,k);
%    s_eci_k = Sun_ECI_ds(:,k);
%    if k>1
%        [R_eci_body_est(:,:,k) P_k]= est_ekf(R_eci_body_est(:,:,k-1), drift_bias, P_k, G_rate_est(:,k), B_sens_est(:,k), I_sun_sens_ds(:,k), norm_sun_sens, b_eci_k, s_eci_k, timestep);
%    else	
%	R_eci_body_est(:,:,k) = est_svd([B_sens_est(:,k),est_sun_sens(I_sun_sens_ds(:,k), norm_sun_sens)],[b_eci_k,s_eci_k]);
%    end

    %SVD Least Squares
    %Assumes start in sun
    if k>1
    	w_ss = cart2ss(G_rate_est(:,k));
        diff_x_est = w_ss*R_test(:,:,k-1);
        R_test(:,:,k) = R_eci_body_est(:,:,k-1) + diff_x_est*timestep;
    end

    b_eci_k = B_ECI_ds(:,k);
    s_eci_k = Sun_ECI_ds(:,k);
    ss_cnt = 0;
    for ss = 1:length(I_sun_sens_ds(:,1))
        if (I_sun_sens_ds(ss,k))
            ss_cnt = ss_cnt + 1;
        end
    end
    if (ss_cnt >= 3)
        R_eci_body_est(:,:,k) = est_svd([B_sens_est(:,k),est_sun_sens(I_sun_sens_ds(:,k), norm_sun_sens)],[b_eci_k,s_eci_k]);
        %R_eci_body_est(:,:,k) = est_quest([B_sens_est(:,k),est_sun_sens(I_sun_sens_ds(:,k), norm_sun_sens)],[b_eci_k,s_eci_k]);
    else
        x_w_int = R_test(:,:,k)'*[1;0;0];
        R_eci_body_est(:,:,k) = est_svd([B_sens_est(:,k),x_w_int],[b_eci_k,[1;0;0]]);
        %R_eci_body_est(:,:,k) = est_quest([B_sens_est(:,k),x_w_int],[b_eci_k,[1;0;0]]);
    end

    fprintf(1,'\b\b\b%02d%%',floor((k/sens_length)*100));
end
fprintf(1,'\n');

Sat_ECI_xhat_est = (reshape(R_eci_body_est(:,1,:),3,sens_length));
Sat_ECI_yhat_est = (reshape(R_eci_body_est(:,2,:),3,sens_length));
Sat_ECI_zhat_est = (reshape(R_eci_body_est(:,3,:),3,sens_length));

%% Metrics
fprintf('Metrics\n');
for i=1:length(Sat_ECI_zhat(1,:))
    find_time_ind = find(Time_rs <= Orbit.Time(i));
    zhat_est = Sat_ECI_zhat_est(:,find_time_ind(length(find_time_ind)));
    xhat_est = Sat_ECI_xhat_est(:,find_time_ind(length(find_time_ind)));
    ram_err(:,i) = acosd(dot(zhat_est,Sat_ECI_zhat(:,i))./(norm(zhat_est)*norm(Sat_ECI_zhat(:,i))));
    sun_err(:,i) = acosd(dot(xhat_est,Sat_ECI_xhat(:,i))./(norm(xhat_est)*norm(Sat_ECI_xhat(:,i))));
end

mean_ram_err = mean(ram_err);
var_ram_err = var(ram_err);
sigma_ram_err = sqrt(var_ram_err);
mean_sun_err = mean(sun_err);
var_sun_err = var(sun_err);
sigma_sun_err = sqrt(var_sun_err);

deg_v = 0:.01:max(ram_err);
ram_err_1_sigma = 0;
ram_err_3_sigma = 0;
for i = 1:length(deg_v)
    ram_err_cdf(:,i) = nnz(ram_err(:,:) < deg_v(i))/length(ram_err);
    if ram_err_cdf(:,i) >= .6827 && ram_err_1_sigma == 0
        ram_err_1_sigma = deg_v(i);
    end
    if ram_err_cdf(:,i) >= .9973 && ram_err_3_sigma == 0
        ram_err_3_sigma = deg_v(i);
    end
end

fprintf('Mean Ram Error:          %2.2f Deg     Sigma: %2.2f Deg\n',mean_ram_err,ram_err_1_sigma);
fprintf('Mean Sun Pointing Error: %2.2f Deg     \n',mean_sun_err);

%% Plots
%%Plots
figure;
plot(deg_v,ram_err_cdf);
title('Ram Error Cumulative Histogram');
hold on;
yL = get(gca,'YLim');
line([ram_err_1_sigma,ram_err_1_sigma],yL,'Color','g');
line([ram_err_3_sigma,ram_err_3_sigma],yL,'Color','r');
xlabel('Error (deg)');

%Decimate Orbit Sim vectors for easier plotting
for i=1:length(Orbit.Time)
	if mod(i,10) == 0
		Orbit_Time_dec(:,i/10) = Orbit.Time(i);
		ram_err_dec(:,i/10) = ram_err(i);
		sun_err_dec(:,i/10) = sun_err(i);
                Sat_ECI_xrate_dec(:,i/10) = Sat_ECI_xrate(i);
                Sat_ECI_yrate_dec(:,i/10) = Sat_ECI_yrate(i);
                Sat_ECI_zrate_dec(:,i/10) = Sat_ECI_zrate(i);
		Orbit_Sun_ECI_dec(:,i/10) = Orbit.Sun_ECI(:,i);
		S_true_bf_dec(:,i/10) = S_true_bf(:,i);
		B_true_bf_dec(:,i/10) = B_true_bf(:,i);
		Orbit_B_ECI_dec(:,i/10) = Orbit.B_ECI(:,i);
	end
end

B_true_norm = sqrt(sum((B_true_bf).^2));
S_true_norm = sqrt(sum((S_true_bf).^2));
B_sens_norm = sqrt(sum((B_sens_est).^2));
%S_sens_norm = sqrt(sum((S_sens_est_bf).^2));
%B_err = sqrt(sum((B_true_bf-B_sens_est).^2));
%S_err = sqrt(sum((S_true_bf-S_sens_est_bf).^2));
B_angle_err = acosd(dot(B_true_bf(:,1:length(B_sens_est(1,:))),B_sens_est)./(B_true_norm(:,1:length(B_sens_est(1,:))).*B_sens_norm));
%S_angle_err = acosd(dot(S_true_bf,S_sens_est_bf)./(S_true_norm.*S_sens_norm));

%Plot Error
figure;
subplot(2,1,1);
plot(Orbit_Time_dec,ram_err_dec);
title('Ram Pointing Error');
hold on;
plot(Orbit_Time_dec,2*ones(length(Orbit_Time_dec)),'r');
ylabel('Error (degrees)');
xlabel('Time (s)');
subplot(2,1,2);
plot(Orbit_Time_dec,sun_err_dec);
title('Sun Pointing Error');
hold on;
plot(Orbit_Time_dec,2*ones(length(Orbit_Time_dec)),'r');
ylabel('Error (degrees)');
xlabel('Time (s)');


%%Plot Body Rates
figure;
subplot(2,1,1);
plot(Orbit_Time_dec,(180/pi())*Sat_ECI_xrate_dec,'r');
title('Body Rates');
hold on;
plot(Orbit_Time_dec,(180/pi())*Sat_ECI_yrate_dec,'g');
plot(Orbit_Time_dec,(180/pi())*Sat_ECI_zrate_dec,'b');
ylabel('Deg/s');
xlabel('Time (s)');
hold off;
subplot(2,1,2);
plot(Time_rs,(180/pi())*G_rate_est(1,:),'r');
title('Body Rates from Gyro');
ylabel('Deg/s');
xlabel('Time (s)');
hold on;
plot(Time_rs,(180/pi())*G_rate_est(2,:),'g');
plot(Time_rs,(180/pi())*G_rate_est(3,:),'b');
hold off

%Plot Magnetometer
figure;
subplot(3,1,1);
hold on;
plot(Orbit_Time_dec,B_true_bf_dec(1,:),'r');
plot(Orbit_Time_dec,B_true_bf_dec(2,:),'g');
plot(Orbit_Time_dec,B_true_bf_dec(3,:),'b');
title('Magnetic Field Vector');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(3,1,2);
plot(Time_rs,B_sens_est(1,:),'r');
hold on;
plot(Time_rs,B_sens_est(2,:),'g');
plot(Time_rs,B_sens_est(3,:),'b');
title('Magnetic Field Vector Estimate');
legend('B_x','B_y','B_z');
xlabel('Time');
ylabel('Magnetic Field Strength (Magnetometer)')
hold off;
subplot(3,1,3);
plot(Time_rs, B_angle_err);
title('Magnetic Field Vector Error');
ylabel('Degrees');
xlabel('Time');

%Sun Estimate
figure;
subplot(4,1,1);
hold on;
plot(Orbit_Time_dec,S_true_bf_dec(1,:),'r');
plot(Orbit_Time_dec,S_true_bf_dec(2,:),'g');
plot(Orbit_Time_dec,S_true_bf_dec(3,:),'b');
title('Sun Vector');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(4,1,2);
hold on;
%plot(Time_rs,S_sens_est_bf(1,:),'r');
%plot(Time_rs,S_sens_est_bf(2,:),'g');
%plot(Time_rs,S_sens_est_bf(3,:),'b');
title('Sun Vector Estimate');
legend('BF\_X','BF\_Y','BF\_Z');
hold off;
%subplot(4,1,3);
%plot(Time_rs, S_sens_num_bf);
%title('Number of active Sun sensors');
%subplot(4,1,4);
%plot(Time_rs, S_angle_err);
%title('Sun Vector Error');

%Body Frame
figure;
subplot(2,3,1);
hold on;
plot(Orbit.Time,Sat_ECI_xhat(1,:),'r');
plot(Orbit.Time,Sat_ECI_xhat(2,:),'g');
plot(Orbit.Time,Sat_ECI_xhat(3,:),'b');
title('Satellite x\_hat (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(2,3,2);
hold on;
plot(Orbit.Time,Sat_ECI_yhat(1,:),'r');
plot(Orbit.Time,Sat_ECI_yhat(2,:),'g');
plot(Orbit.Time,Sat_ECI_yhat(3,:),'b');
title('Satellite y\_hat (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(2,3,3);
hold on;
plot(Orbit.Time,Sat_ECI_zhat(1,:),'r');
plot(Orbit.Time,Sat_ECI_zhat(2,:),'g');
plot(Orbit.Time,Sat_ECI_zhat(3,:),'b');
title('Satellite z\_hat (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(2,3,4);
hold on;
plot(Time_rs,Sat_ECI_xhat_est(1,:),'r');
plot(Time_rs,Sat_ECI_xhat_est(2,:),'g');
plot(Time_rs,Sat_ECI_xhat_est(3,:),'b');
title('Satellite x\_hat Estimate (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
axis([0,length(Time_rs)/2,-1,1]);
hold off;
subplot(2,3,5);
hold on;
plot(Time_rs,Sat_ECI_yhat_est(1,:),'r');
plot(Time_rs,Sat_ECI_yhat_est(2,:),'g');
plot(Time_rs,Sat_ECI_yhat_est(3,:),'b');
title('Satellite y\_hat estimate (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
axis([0,length(Time_rs)/2,-1,1]);
hold off;
subplot(2,3,6);
hold on;
plot(Time_rs,Sat_ECI_zhat_est(1,:),'r');
plot(Time_rs,Sat_ECI_zhat_est(2,:),'g');
plot(Time_rs,Sat_ECI_zhat_est(3,:),'b');
title('Satellite z\_hat estimate (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
axis([0,length(Time_rs)/2,-1,1]);
hold off;
;

