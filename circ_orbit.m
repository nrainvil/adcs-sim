%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nicholas Rainville, Jacob Cook
% 10/9/2013
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear all;
close all;

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
OE.i        = deg2rad(98);      % inclination [rad]
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
startTime = datenum('01 May 2015 00:00:00');
stopTime  = datenum('01 May 2015 04:00:00');
timeStep = 10; % [sec]
time = startTime:timeStep/(3600*24):stopTime;

Orbit = GenOrbit(time, TLE);
clear SatNum SatClass IntDsgntr Epoch startTime stopTime timeStep time

%% Create Satellite Frame
Sat_ECI_xhat = zeros(size(Orbit.posECI_km));
Sat_ECI_zhat = zeros(size(Orbit.posECI_km));
Sat_ECI_yhat = zeros(size(Orbit.posECI_km));

Sat_ECI_zhat(1,:) = -1*Orbit.velECI_kmps(1,:)./sqrt(sum(abs(Orbit.velECI_kmps).^2,1));
Sat_ECI_zhat(2,:) = -1*Orbit.velECI_kmps(2,:)./sqrt(sum(abs(Orbit.velECI_kmps).^2,1));
Sat_ECI_zhat(3,:) = -1*Orbit.velECI_kmps(3,:)./sqrt(sum(abs(Orbit.velECI_kmps).^2,1));

for i=1:length(Orbit.velECI_kmps(1,:))
    Sat_ECI_xhat(:,i) = ortho_norm(-Orbit.posECI_km(:,i),Sat_ECI_zhat(:,i));
    Sat_ECI_yhat(:,i) = cross(Sat_ECI_zhat(:,i),Sat_ECI_xhat(:,i));
    B_norm(:,i) = Orbit.B_ECI(:,i)/norm(Orbit.B_ECI(:,i));
end

%% Sensors

%3x Magnetometers
B_Mag_1 = mag_sens(Orbit.B_ECI,Sat_ECI_xhat);
B_Mag_2 = mag_sens(Orbit.B_ECI,Sat_ECI_yhat);
B_Mag_3 = mag_sens(Orbit.B_ECI,Sat_ECI_zhat);

%14x Sun sensors (RAX2 positioning)
norm_sun_sens(1,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 17,-10);
norm_sun_sens(2,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,20);
norm_sun_sens(3,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -17,10);
norm_sun_sens(4,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -162,-10);
norm_sun_sens(5,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 180,20);
norm_sun_sens(6,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 162,-10);
norm_sun_sens(7,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 72,10);
norm_sun_sens(8,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 107,10);
norm_sun_sens(9,:,:)  = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 90,-20);
norm_sun_sens(10,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -107,10);
norm_sun_sens(11,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -72,10);
norm_sun_sens(12,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, -90,-20);
norm_sun_sens(13,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,90);
norm_sun_sens(14,:,:) = rot_azel(Sat_ECI_xhat, Sat_ECI_zhat, 0,-90);

[I_sun_sens, I_sun_sens_alpha]   = sun_sens(Orbit.Sun_ECI,norm_sun_sens); 

%Estimate sun vector solution based on sensor
[S_sens_est, S_sens_num] = est_sun_sens(I_sun_sens, norm_sun_sens);

%% Generate Cosine Matrix
% for i=1:length(B_Mag_1)
%     b_b = [B_Mag_1(:,i); B_Mag_2(:,i); B_Mag_3(:,i)];
%     R_bi = triad(Orbit.Sun_ECI(:,i),Orbit.Sun_ECI(:,i),b_b,Orbit.B_ECI(:,i));
% end


%% Plots
figure;
subplot(3,1,1);
hold on;
plot(Sat_ECI_xhat(1,:),'r');
plot(Sat_ECI_xhat(2,:),'g');
plot(Sat_ECI_xhat(3,:),'b');
title('Satellite x\_hat (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(3,1,2);
hold on;
plot(Sat_ECI_yhat(1,:),'r');
plot(Sat_ECI_yhat(2,:),'g');
plot(Sat_ECI_yhat(3,:),'b');
title('Satellite y\_hat (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(3,1,3);
hold on;
plot(Sat_ECI_zhat(1,:),'r');
plot(Sat_ECI_zhat(2,:),'g');
plot(Sat_ECI_zhat(3,:),'b');
title('Satellite z\_hat (ECI)')
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;

%Plot Magnetometer
figure;
subplot(2,1,1);
hold on;
plot(Orbit.B_ECI(1,:),'r');
plot(Orbit.B_ECI(2,:),'g');
plot(Orbit.B_ECI(3,:),'b');
title('Magnetic Field Vector');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(2,1,2);
plot(Orbit.Time,B_Mag_1,'r');
hold on;
plot(Orbit.Time,B_Mag_2,'g');
plot(Orbit.Time,B_Mag_3,'b');
legend('B_x','B_y','B_z');
xlabel('Time');
ylabel('Magnetic Field Strength (Magnetometer)')
hold off;

%Sun Estimate
figure;
subplot(3,1,1);
hold on;
plot(Orbit.Time,Orbit.Sun_ECI(1,:),'r');
plot(Orbit.Time,Orbit.Sun_ECI(2,:),'g');
plot(Orbit.Time,Orbit.Sun_ECI(3,:),'b');
title('Sun Vector');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(3,1,2);
hold on;
plot(Orbit.Time,S_sens_est(1,:),'r');
plot(Orbit.Time,S_sens_est(2,:),'g');
plot(Orbit.Time,S_sens_est(3,:),'b');
title('Sun Vector Estimate');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(3,1,3);
plot(Orbit.Time, S_sens_num);
title('Number of active Sun sensors');

%Plot Sun Sensor
figure;
subplot(2,3,1);
plot(Orbit.Time,I_sun_sens(1,:),'r');
hold on;
plot(Orbit.Time,I_sun_sens(2,:),'g');
plot(Orbit.Time,I_sun_sens(3,:),'c');
hold off;
title('+x');
subplot(2,3,4);
plot(Orbit.Time,I_sun_sens(4,:),'r');
hold on;
plot(Orbit.Time,I_sun_sens(5,:),'g');
plot(Orbit.Time,I_sun_sens(6,:),'c');
hold off;
title('-x');
subplot(2,3,2);
plot(Orbit.Time,I_sun_sens(7,:),'r');
hold on;
plot(Orbit.Time,I_sun_sens(8,:),'g');
plot(Orbit.Time,I_sun_sens(9,:),'c');
hold off;
title('+y');
subplot(2,3,5);
plot(Orbit.Time,I_sun_sens(10,:),'r');
hold on;
plot(Orbit.Time,I_sun_sens(11,:),'g');
plot(Orbit.Time,I_sun_sens(12,:),'c');
hold off;
title('-y');
subplot(2,3,3);
plot(Orbit.Time,I_sun_sens(13,:),'r');
hold off;
title('+z');
subplot(2,3,6);
plot(Orbit.Time,I_sun_sens(14,:),'r');
hold off;
title('-z');



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
