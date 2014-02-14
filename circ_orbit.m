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

Sat_ECI_xhat(1,:) = -1*Orbit.posECI_km(1,:)./sqrt(sum(abs(Orbit.posECI_km).^2,1));
Sat_ECI_xhat(2,:) = -1*Orbit.posECI_km(2,:)./sqrt(sum(abs(Orbit.posECI_km).^2,1));
Sat_ECI_xhat(3,:) = -1*Orbit.posECI_km(3,:)./sqrt(sum(abs(Orbit.posECI_km).^2,1));

for i=1:length(Orbit.velECI_kmps(1,:))
    Sat_ECI_zhat(:,i) = ortho_norm(Orbit.velECI_kmps(:,i),Sat_ECI_xhat(:,i));
    Sat_ECI_yhat(:,i) = cross(Sat_ECI_zhat(:,i),Sat_ECI_xhat(:,i));
    B_norm(:,i) = Orbit.B_ECI(:,i)/norm(Orbit.B_ECI(:,i));
end

%% Sensors
I_0 = 1/norm(Orbit.Sun_ECI(:,1));
Sens_Sun_1_normal = -1*Sat_ECI_xhat;
for i = 1:length(Orbit.Sun_ECI(1,:))
    dot_alpha_1 = dot(Orbit.Sun_ECI(:,i),Sens_Sun_1_normal(:,i))/norm(Orbit.Sun_ECI(:,i));
    Sens_Sun_1_alpha(1,i) = acosd(dot_alpha_1)*dot_alpha_1/abs(dot_alpha_1);
    
    Sun_SAT(1,i) = dot(Orbit.Sun_ECI(:,i),Sat_ECI_xhat(:,i));
    Sun_SAT(2,i) = dot(Orbit.Sun_ECI(:,i),Sat_ECI_yhat(:,i));
    Sun_SAT(3,i) = dot(Orbit.Sun_ECI(:,i),Sat_ECI_zhat(:,i));
    if Sens_Sun_1_alpha(1,i) > -35 && Sens_Sun_1_alpha(1,i) < 35
        I_Sun_1(1,i) = I_0*cosd(Sens_Sun_1_alpha(1,i));
    else
        I_Sun_1(1,i) = 0;
    end
end

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
figure;
title('T0');
hold on;
plot3([0 1],[0 0],[0 0],'r',[0 0],[0 1],[0 0],'g',[0 0],[0 0],[0 1],'b','linewidth',2);
plot3([0 Sat_ECI_xhat(1,1)],[0 Sat_ECI_xhat(2,1)],[0 Sat_ECI_xhat(3,1)],'r');
plot3([0 Sat_ECI_yhat(1,1)],[0 Sat_ECI_yhat(2,1)],[0 Sat_ECI_yhat(3,1)],'g');
plot3([0 Sat_ECI_zhat(1,1)],[0 Sat_ECI_zhat(2,1)],[0 Sat_ECI_zhat(3,1)],'r');
plot3([0 B_norm(1,1)],[0 B_norm(2,1)],[0 B_norm(3,1)],'c');
plot3([0 Orbit.Sun_ECI(1,1)],[0 Orbit.Sun_ECI(2,1)],[0 Orbit.Sun_ECI(3,1)],'k');
xlabel('X');ylabel('Y');zlabel('Z');
set(gca,'DataAspectRatio',[1 1 1]);
view([1,1,1]);
hold off;

%Plot Sun Sensor
figure;
subplot(3,1,1);
plot(I_Sun_1);
title('I @ 0 deg');
subplot(3,1,2);
plot(Sens_Sun_1_alpha);
title('alpha angle');
subplot(3,1,3);
title('Sun Vector Sat Body Frame');
hold on;
plot(Sun_SAT(1,:),'r');
plot(Sun_SAT(2,:),'g');
plot(Sun_SAT(3,:),'b');
legend('Sat Body X\_hat','Sat Body Y\_hat', 'Sat_Body Z\_hat');
hold off;

%% Animation
fig_anim = figure('Position',[100,100,900,900]);
for k = 1:length(Sat_ECI_xhat(1,:))
    plot3(Sat_ECI_xhat(1,k),Sat_ECI_xhat(2,k),Sat_ECI_xhat(3,k),'r','marker','o'); 
    hold on;
    plot3([Sat_ECI_xhat(1,k), -Sat_ECI_zhat(1,k)+Sat_ECI_xhat(1,k)],[Sat_ECI_xhat(2,k), -Sat_ECI_zhat(2,k)+Sat_ECI_xhat(2,k)],[Sat_ECI_xhat(3,k), -Sat_ECI_zhat(3,k)+Sat_ECI_xhat(3,k)],'b','linewidth',2);
    plot3([0 1/2],[0 0],[0 0],'r',[0 0],[0 1/2],[0 0],'g',[0 0],[0 0],[0 1/2],'b','linewidth',2);
    plot3([Sat_ECI_xhat(1,k), -B_norm(1,k)+Sat_ECI_xhat(1,k)],[Sat_ECI_xhat(2,k), -B_norm(2,k)+Sat_ECI_xhat(2,k)],[Sat_ECI_xhat(3,k), -B_norm(3,k)+Sat_ECI_xhat(3,k)],'c','linewidth',2);
    plot3([Sat_ECI_xhat(1,k), -Orbit.Sun_ECI(1,k)+Sat_ECI_xhat(1,k)],[Sat_ECI_xhat(2,k), -Orbit.Sun_ECI(2,k)+Sat_ECI_xhat(2,k)],[Sat_ECI_xhat(3,k) -Orbit.Sun_ECI(3,k)+Sat_ECI_xhat(3,k)],'y','linewidth',2);
    xlabel('X');ylabel('Y');zlabel('Z');
    axis([-1,1,-1,1,-1,1]);
    view([1,0,0]);
    hold off;
    M(k) = getframe; 
end
