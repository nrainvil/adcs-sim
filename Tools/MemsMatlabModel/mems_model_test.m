clear all; close all; clc

% Load real gyro measurements
load('gyro_data_test.mat');
real_gyro_data = xg(:,1)';

%% MEMS PARAMETERS
%--- Set MEMS parameters --------------------------------------------------
samp_freq          = 50;          % gyro sampling frequency [Hz]
% samp_freq          = 95;
sensitivity        = 0.010;     % 8.75 mdps/digit
% sensitivity        = 0.00875;
gyro_bias          = 0;           % Static bias (dps)
rate_noise_density = 0.014;        % [dps/sqrt(Hz)] (fft)
% rate_noise_density = 0.03;
upper_saturation   = 2^15;        % 16 bits resolution
lower_saturation   = -2^15;   
wn                 = 1000*2*pi;   % MEMS dynamics [rad/s]
w                  = 12.5*2*pi;   % LPF cut off frequency [rad/s]
zeta               = sqrt(2)/2;   % Damping estimate
z                  = zeta;

%--- Plot raw real gyro data
time = 0 : 1/samp_freq : (length(real_gyro_data)-1)/samp_freq;
figure
plot(time,real_gyro_data)
title('Real gyro measurement')
xlabel('Time (sec)')
ylabel('Volts')

%% MEMS CALIBRATION (used only when hardware is available to test)
%--- Find bias ------------------------------------------------------------
v2dps     = sensitivity;
gyro_bias = mean(real_gyro_data) * v2dps; 
% gyro_bias = -0.6851;  % dps
real_gyro_data_noBias = real_gyro_data - mean(real_gyro_data);

%--- Calibrate Rate Noise Density -----------------------------------------
[tau, AVAR] = myAllan(time, real_gyro_data_noBias*v2dps);
[index value] = find(tau == 1)
rate_noise_density = AVAR(samp_freq);
% rate_noise_density = 0.0375;

%% MEMS model
real_rate = [time', zeros(length(time),1)];      % Either real world body rate or simulated (0 = stationary)
sim_tim   = time(end);   % Length of the simulation [sec]

sim_gyro      = sim('MEMS_model','StopTime','sim_tim','SaveOutput','on','FixedStep','1/samp_freq');
sim_gyro_data = sim_gyro.get('sensed_angular_rate');

%% PLOT RESULTS
% Time domain sim gyro data
figure
plot(time, real_gyro_data, 'r')
hold all
plot(time, sim_gyro_data, 'b')
ylabel('Sensed Angular Rate [Volts]')
xlabel('Time [sec]')
title('Real and Simulated Gyro Measurements')
legend('Real', 'Sim')

% Frequency content
sim_gyro_noBias = sim_gyro_data - mean(sim_gyro_data);  % Remove bias from the simulated gyro data
figure
[f, amp] = myfft(real_gyro_data_noBias * v2dps, 1/samp_freq);
hold all
[f, amp] = myfft(sim_gyro_noBias * v2dps, 1/samp_freq);
legend('Real', 'Sim')

















