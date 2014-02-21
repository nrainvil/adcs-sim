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

for i = 1:length(Sat_ECI_xhat)
    norm_sun_sens_1(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 17,-10);
    norm_sun_sens_2(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 0,20);
    norm_sun_sens_3(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), -17,10);
    norm_sun_sens_4(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), -162,-10);
    norm_sun_sens_5(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 180,20);
    norm_sun_sens_6(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 162,-10);
    norm_sun_sens_7(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 72,10);
    norm_sun_sens_8(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 107,10);
    norm_sun_sens_9(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 90,-20);
    norm_sun_sens_10(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), -107,10);
    norm_sun_sens_11(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), -72,10);
    norm_sun_sens_12(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), -90,-20);
    norm_sun_sens_13(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 0,90);
    norm_sun_sens_14(:,i) = rot_azel(Sat_ECI_xhat(:,i), Sat_ECI_zhat(:,i), 0,-90);
end

[I_Sun_1, I_Sun_alpha_1] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_1); 
[I_Sun_2, I_Sun_alpha_2] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_2);
[I_Sun_3, I_Sun_alpha_3] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_3);
[I_Sun_4, I_Sun_alpha_4] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_4); 
[I_Sun_5, I_Sun_alpha_5] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_5);
[I_Sun_6, I_Sun_alpha_6] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_6);
[I_Sun_7, I_Sun_alpha_7] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_7); 
[I_Sun_8, I_Sun_alpha_8] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_8);
[I_Sun_9, I_Sun_alpha_9] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_9);
[I_Sun_10, I_Sun_alpha_10] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_10); 
[I_Sun_11, I_Sun_alpha_11] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_11);
[I_Sun_12, I_Sun_alpha_12] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_12);
[I_Sun_13, I_Sun_alpha_13] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_13);
[I_Sun_14, I_Sun_alpha_14] = sun_sens(Orbit.Sun_ECI,norm_sun_sens_14);

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

%Sun and Mag
figure;
subplot(2,1,1);
hold on;
plot(Orbit.Sun_ECI(1,:),'r');
plot(Orbit.Sun_ECI(2,:),'g');
plot(Orbit.Sun_ECI(3,:),'b');
title('Sun Vector');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;
subplot(2,1,2);
hold on;
plot(Orbit.B_ECI(1,:),'r');
plot(Orbit.B_ECI(2,:),'g');
plot(Orbit.B_ECI(3,:),'b');
title('Magnetic Field Vector');
legend('ECI\_X','ECI\_Y','ECI\_Z');
hold off;

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

%Plot Sun Sensor
figure;
subplot(2,3,1);
plot(Orbit.Time,I_Sun_1,'r');
hold on;
plot(Orbit.Time,I_Sun_2,'g');
plot(Orbit.Time,I_Sun_3,'c');
hold off;
title('+x');
subplot(2,3,4);
plot(Orbit.Time,I_Sun_4,'r');
hold on;
plot(Orbit.Time,I_Sun_5,'g');
plot(Orbit.Time,I_Sun_6,'c');
hold off;
title('-x');
subplot(2,3,2);
plot(Orbit.Time,I_Sun_7,'r');
hold on;
plot(Orbit.Time,I_Sun_8,'g');
plot(Orbit.Time,I_Sun_9,'c');
hold off;
title('+y');
subplot(2,3,5);
plot(Orbit.Time,I_Sun_10,'r');
hold on;
plot(Orbit.Time,I_Sun_11,'g');
plot(Orbit.Time,I_Sun_12,'c');
hold off;
title('-y');
subplot(2,3,3);
plot(Orbit.Time,I_Sun_13,'r');
hold off;
title('+z');
subplot(2,3,6);
plot(Orbit.Time,I_Sun_13,'r');
hold off;
title('-z');


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
