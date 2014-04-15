clear all;
close all;


t = 0:.1:100;
real_gyro_data = sin((2*pi/10)*(t.^2/100));
samp_freq = 1; %Hz

R = 1000;
C = .001;

real_rate = [t',real_gyro_data'];
sim_gyro = sim('RC');
%sim_gyro_data = sim_gyro.get('sensed_angular_rate');

figure;
plot(t,real_gyro_data,'r');
hold on;
plot(sim_gyro,sensed_angular_rate,'b');
