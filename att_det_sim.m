%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nicholas Rainville, Jacob Cook
% 10/9/2013
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all;
close all;
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

%% Set Classical Orbital Elements
% Note: We start at closest point to south pole to line up with our K
OE.a        = 370 + R_EARTH;    % semi-major axis [km]
OE.e        = 0;                % eccentricity
OE.i        = deg2rad(90);      % inclination [rad]
OE.omega    = deg2rad(0);       % argument of perigee [rad]
OE.RAAN     = deg2rad(30);       % right ascension of the ascending node [rad]
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
stopTime  = datenum('01 May 2015 01:32:00');
timeStep = .2; % [sec] (5 Hz)
time = startTime:timeStep/(3600*24):stopTime;

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

%% ECI Body Rates
for i=2:length(Sat_ECI_xhat)
    R_m = reshape(R_eci_body(:,:,i-1),3,3);
    R_p = reshape(R_eci_body(:,:,i),3,3);
    dr_dt = ((R_p-R_m)/timeStep);
    w_i = dr_dt*R_p';
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
[B_mag_sens, I_sun_sens, G_rate_sens, norm_sun_sens] = sc_sens_sim(R_eci_body, R_dot_eci_body, Orbit.Sun_ECI,Orbit.B_ECI);

%Discrete
%[B_mag_sens_rs, I_sun_sens_rs, G_rate_sens_rs, Time_rs]
fs_sens = 1/20; %Hz
Sun_ECI_ts = timeseries(Orbit.Sun_ECI,Orbit.Time,'name','I_eci');
I_sun_sens_ts = timeseries(I_sun_sens,Orbit.Time,'name','I_bf');
B_ECI_ts = timeseries(Orbit.B_ECI,Orbit.Time,'name','B_eci');
B_mag_sens_ts = timeseries(B_mag_sens,Orbit.Time,'name','B_bf');
G_rate_sens_ts = timeseries(G_rate_sens,Orbit.Time,'name','G_bf');

Time_rs = Orbit.Time(1):(1/fs_sens):(Orbit.Time(length(Orbit.Time)));
Sun_ECI_rs = resample(Sun_ECI_ts,Time_rs);
I_sun_sens_rs = resample(I_sun_sens_ts,Time_rs);
B_ECI_rs = resample(B_ECI_ts,Time_rs);
B_mag_sens_rs = resample(B_mag_sens_ts,Time_rs);
G_rate_sens_rs = resample(G_rate_sens_ts,Time_rs);
sens_length = length(Time_rs);

%% Attitude Estimate
%Estimate sun vector
I_sun_sens_ds = reshape(I_sun_sens_rs.data,length(I_sun_sens(:,1)),sens_length);
[S_sens_est_bf, S_sens_num_bf] = est_sun_sens(I_sun_sens_ds, norm_sun_sens);

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
for k=1:sens_length    
    if k>1
        w_ss = [0,-G_rate_est(3,k),G_rate_est(2,k); ...
                G_rate_est(3,k), 0, -G_rate_est(1,k); ...
                -G_rate_est(2,k), G_rate_est(1,k), 0;];
        diff_x_est = w_ss*R_test(:,:,k-1);
        R_test(:,:,k) = R_eci_body_est(:,:,k-1) + diff_x_est*timestep;
    end
    b_k = [B_sens_est(:,k),S_sens_est_bf(:,k),[1;0;0],[0;0;1]];
    eci_k = [B_ECI_ds(:,k),Sun_ECI_ds(:,k),R_test(:,1,k),R_test(:,3,k)];
    R_eci_body_est(:,:,k) = est_svd(b_k, eci_k);
end

Sat_ECI_xhat_est = (reshape(R_eci_body_est(:,1,:),3,sens_length));
Sat_ECI_yhat_est = (reshape(R_eci_body_est(:,2,:),3,sens_length));
Sat_ECI_zhat_est = (reshape(R_eci_body_est(:,3,:),3,sens_length));

%% Metrics
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

%figure;
%subplot(2,1,1);
%plot(Orbit.Time,ram_err);
%title('Ram Pointing Error');
%hold on;
%plot(Orbit.Time,2*ones(length(Orbit.Time)),'r');
%ylabel('Error (degrees)');
%xlabel('Time (s)');
%subplot(2,1,2);
%plot(Orbit.Time,sun_err);
%title('Sun Pointing Error');
%hold on;
%plot(Orbit.Time,2*ones(length(Orbit.Time)),'r');
%ylabel('Error (degrees)');
%xlabel('Time (s)');

% figure;hist(ram_err,50);hold on;
% yL = get(gca,'YLim');
% line([mean_ram_err,mean_ram_err],yL,'Color','g');
% line([mean_ram_err+sigma_ram_err,mean_ram_err+sigma_ram_err],yL,'Color','r');
% xlabel('Error (degrees)');
% ylabel('Count');

%figure;
%subplot(2,3,1);
%hold on;
%plot(Orbit.Time,Sat_ECI_xhat(1,:),'r');
%plot(Orbit.Time,Sat_ECI_xhat(2,:),'g');
%plot(Orbit.Time,Sat_ECI_xhat(3,:),'b');
%title('Satellite x\_hat (ECI)')
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%subplot(2,3,2);
%hold on;
%plot(Orbit.Time,Sat_ECI_yhat(1,:),'r');
%plot(Orbit.Time,Sat_ECI_yhat(2,:),'g');
%plot(Orbit.Time,Sat_ECI_yhat(3,:),'b');
%title('Satellite y\_hat (ECI)')
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%subplot(2,3,3);
%hold on;
%plot(Orbit.Time,Sat_ECI_zhat(1,:),'r');
%plot(Orbit.Time,Sat_ECI_zhat(2,:),'g');
%plot(Orbit.Time,Sat_ECI_zhat(3,:),'b');
%title('Satellite z\_hat (ECI)')
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%subplot(2,3,4);
%hold on;
%plot(Time_rs,Sat_ECI_xhat_est(1,:),'r');
%plot(Time_rs,Sat_ECI_xhat_est(2,:),'g');
%plot(Time_rs,Sat_ECI_xhat_est(3,:),'b');
%title('Satellite x\_hat Estimate (ECI)')
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%subplot(2,3,5);
%hold on;
%plot(Time_rs,Sat_ECI_yhat_est(1,:),'r');
%plot(Time_rs,Sat_ECI_yhat_est(2,:),'g');
%plot(Time_rs,Sat_ECI_yhat_est(3,:),'b');
%title('Satellite y\_hat estimate (ECI)')
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%subplot(2,3,6);
%hold on;
%plot(Time_rs,Sat_ECI_zhat_est(1,:),'r');
%plot(Time_rs,Sat_ECI_zhat_est(2,:),'g');
%plot(Time_rs,Sat_ECI_zhat_est(3,:),'b');
%title('Satellite z\_hat estimate (ECI)')
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%
%%Plot Body Rates
%figure;
%subplot(2,1,1);
%plot(Orbit.Time,(180/pi())*Sat_ECI_xrate,'r');
%title('Body Rates');
%hold on;
%plot(Orbit.Time,(180/pi())*Sat_ECI_yrate,'g');
%plot(Orbit.Time,(180/pi())*Sat_ECI_zrate,'b');
%ylabel('Deg/s');
%xlabel('Time (s)');
%hold off;
%subplot(2,1,2);
%plot(Time_rs,(180/pi())*G_rate_est(1,:),'r');
%title('Body Rates from Gyro');
%ylabel('Deg/s');
%xlabel('Time (s)');
%hold on;
%plot(Time_rs,(180/pi())*G_rate_est(2,:),'g');
%plot(Time_rs,(180/pi())*G_rate_est(3,:),'b');
%hold off;
%
%%Plot Rate Gyro 
%figure;
%subplot(2,1,1);
%hold on;
%plot(Orbit.Time,reshape(R_eci_body(1,3,:),1,length(R_eci_body(1,3,:))),'r');
%plot(Orbit.Time,reshape(R_eci_body(2,3,:),1,length(R_eci_body(2,3,:))),'g');
%plot(Orbit.Time,reshape(R_eci_body(3,3,:),1,length(R_eci_body(3,3,:))),'b');
%title('Satellite x\_hat (ECI)');
%axis tight;
%hold off;
%subplot(2,1,2);
%hold on;
%plot(Time_rs,reshape(R_test(1,3,:),1,length(R_test(1,1,:))),'r');
%plot(Time_rs,reshape(R_test(2,3,:),1,length(R_test(2,1,:))),'g');
%plot(Time_rs,reshape(R_test(3,3,:),1,length(R_test(3,1,:))),'b');
%title('ATTITUDE FROM BODY RATE EST');
%axis tight;
%hold off;
%Plot Magnetometer
%figure;
%subplot(2,1,1);
%hold on;
%plot(Orbit.Time,Orbit.B_ECI(1,:),'r');
%plot(Orbit.Time,Orbit.B_ECI(2,:),'g');
%plot(Orbit.Time,Orbit.B_ECI(3,:),'b');
%title('Magnetic Field Vector');
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%subplot(2,1,2);
%plot(Time_rs,B_sens_est(1,:),'r');
%hold on;
%plot(Time_rs,B_sens_est(2,:),'g');
%plot(Time_rs,B_sens_est(3,:),'b');
%title('Magnetic Field Vector Estimate');
%legend('B_x','B_y','B_z');
%xlabel('Time');
%ylabel('Magnetic Field Strength (Magnetometer)')
%hold off;

%Sun Estimate
%figure;
%subplot(3,1,1);
%hold on;
%plot(Orbit.Time,Orbit.Sun_ECI(1,:),'r');
%plot(Orbit.Time,Orbit.Sun_ECI(2,:),'g');
%plot(Orbit.Time,Orbit.Sun_ECI(3,:),'b');
%title('Sun Vector');
%legend('ECI\_X','ECI\_Y','ECI\_Z');
%hold off;
%% subplot(4,1,2);
%% hold on;
%% plot(Orbit.Time,S_sens_est(1,:),'r');
%% plot(Orbit.Time,S_sens_est(2,:),'g');
%% plot(Orbit.Time,S_sens_est(3,:),'b');
%% title('Sun Vector Estimate');
%% legend('ECI\_X','ECI\_Y','ECI\_Z');
%% hold off;
%subplot(3,1,2);
%hold on;
%plot(Time_rs,S_sens_est_bf(1,:),'r');
%plot(Time_rs,S_sens_est_bf(2,:),'g');
%plot(Time_rs,S_sens_est_bf(3,:),'b');
%title('Sun Vector Estimate');
%legend('BF\_X','BF\_Y','BF\_Z');
%hold off;
%subplot(3,1,3);
%plot(Time_rs, S_sens_num_bf);
%title('Number of active Sun sensors');

%Plot Sun Sensor
%figure;
%subplot(2,3,1);
%plot(Time_rs,I_sun_sens_ds(1,:),'r');
%hold on;
%plot(Time_rs,I_sun_sens_ds(2,:),'g');
%plot(Time_rs,I_sun_sens_ds(3,:),'c');
%ylabel('V');
%xlabel('Time (s)');
%hold off;
%title('+x');
%subplot(2,3,4);
%plot(Time_rs,I_sun_sens_ds(4,:),'r');
%hold on;
%plot(Time_rs,I_sun_sens_ds(5,:),'g');
%plot(Time_rs,I_sun_sens_ds(6,:),'c');
%ylabel('V');
%xlabel('Time (s)');
%hold off;
%title('-x');
%subplot(2,3,2);
%plot(Time_rs,I_sun_sens_ds(7,:),'r');
%hold on;
%plot(Time_rs,I_sun_sens_ds(8,:),'g');
%plot(Time_rs,I_sun_sens_ds(9,:),'c');
%hold off;
%ylabel('V');
%xlabel('Time (s)');
%title('+y');
%subplot(2,3,5);
%plot(Time_rs,I_sun_sens_ds(10,:),'r');
%hold on;
%plot(Time_rs,I_sun_sens_ds(11,:),'g');
%plot(Time_rs,I_sun_sens_ds(12,:),'c');
%ylabel('V');
%xlabel('Time (s)');
%hold off;
%title('-y');
%subplot(2,3,3);
%plot(Time_rs,I_sun_sens_ds(13,:),'r');
%ylabel('V');
%xlabel('Time (s)');
%hold off;
%title('+z');
%subplot(2,3,6);
%plot(Time_rs,I_sun_sens_ds(14,:),'r');
%ylabel('V');
%xlabel('Time (s)');
%hold off;
%title('-z');



%Plot Orbit
% figure;
% title('X\_hat');
% hold on;
% plot3(Sat_ECI_xhat(1,:),Sat_ECI_xhat(2,:),Sat_ECI_xhat(3,:),'r');
% plot3([0 Sat_ECI_yhat(1,1)],[0 Sat_ECI_yhat(2,1)],[0 Sat_ECI_yhat(3,1)],'g');
% plot3(Sat_ECI_zhat(1,:),Sat_ECI_zhat(2,:),Sat_ECI_zhat(3,:),'b');
% plot3([0 1],[0 0],[0 0],'r',[0 0],[0 1],[0 0],'g',[0 0],[0 0],[0 1],'b','linewidth',2);
% plot3(B_norm(1,:),B_norm(2,:),B_norm(3,:),'c');
% plot3([0 Orbit.Sun_ECI(1,1)],[0 Orbit.Sun_ECI(2,1)],[0 Orbit.Sun_ECI(3,1)],'k');
% xlabel('X');ylabel('Y');zlabel('Z');
% set(gca,'DataAspectRatio',[1 1 1]);
% view([1,1,1]);
% hold off;


%%
% figure;
% title('T0');
% hold on;
% plot3([0 1],[0 0],[0 0],'r',[0 0],[0 1],[0 0],'g',[0 0],[0 0],[0 1],'b','linewidth',2);
% plot3([0 Sat_ECI_xhat(1,1)],[0 Sat_ECI_xhat(2,1)],[0 Sat_ECI_xhat(3,1)],'r');
% plot3([0 Sat_ECI_yhat(1,1)],[0 Sat_ECI_yhat(2,1)],[0 Sat_ECI_yhat(3,1)],'g');
% plot3([0 Sat_ECI_zhat(1,1)],[0 Sat_ECI_zhat(2,1)],[0 Sat_ECI_zhat(3,1)],'r');
% plot3([0 B_norm(1,1)],[0 B_norm(2,1)],[0 B_norm(3,1)],'c');
% plot3([0 Orbit.Sun_ECI(1,1)],[0 Orbit.Sun_ECI(2,1)],[0 Orbit.Sun_ECI(3,1)],'k');
% xlabel('X');ylabel('Y');zlabel('Z');
% set(gca,'DataAspectRatio',[1 1 1]);
% view([1,1,1]);
% hold off;



%% Animation
% fig_anim = figure('Position',[100,100,900,900]);
% for k = 1:length(Sat_ECI_xhat(1,:))
%     plot3(Sat_ECI_xhat(1,k),Sat_ECI_xhat(2,k),Sat_ECI_xhat(3,k),'r','marker','o'); 
%     hold on;
%     plot3([Sat_ECI_xhat(1,k), -Sat_ECI_zhat(1,k)+Sat_ECI_xhat(1,k)],[Sat_ECI_xhat(2,k), -Sat_ECI_zhat(2,k)+Sat_ECI_xhat(2,k)],[Sat_ECI_xhat(3,k), -Sat_ECI_zhat(3,k)+Sat_ECI_xhat(3,k)],'b','linewidth',2);
%     plot3([0 1/2],[0 0],[0 0],'r',[0 0],[0 1/2],[0 0],'g',[0 0],[0 0],[0 1/2],'b','linewidth',2);
%     plot3([Sat_ECI_xhat(1,k), -B_norm(1,k)+Sat_ECI_xhat(1,k)],[Sat_ECI_xhat(2,k), -B_norm(2,k)+Sat_ECI_xhat(2,k)],[Sat_ECI_xhat(3,k), -B_norm(3,k)+Sat_ECI_xhat(3,k)],'c','linewidth',2);
%     plot3([Sat_ECI_xhat(1,k), -Orbit.Sun_ECI(1,k)+Sat_ECI_xhat(1,k)],[Sat_ECI_xhat(2,k), -Orbit.Sun_ECI(2,k)+Sat_ECI_xhat(2,k)],[Sat_ECI_xhat(3,k) -Orbit.Sun_ECI(3,k)+Sat_ECI_xhat(3,k)],'y','linewidth',2);
%     xlabel('X');ylabel('Y');zlabel('Z');
%     axis([-1,1,-1,1,-1,1]);
%     view([1,0,0]);
%     hold off;
%     M(k) = getframe; 
% end
