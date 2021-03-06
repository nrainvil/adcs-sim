function [ G_rate ] = rate_sens(true_w)
%RATE_SENS

% %--- Set MEMS parameters --------------------------------------------------
% samp_freq          = 1/step;          % gyro sampling frequency [Hz]
% sensitivity        = 0.010;     % 8.75 mdps/digit
gyro_bias          = 0.0;           % Static bias (dps)
rate_noise_density = .05; %.014; %0.014;        % [dps/sqrt(Hz)] (fft) (Epson .004) (ADIS .05)
% upper_saturation   = 2^15;        % 16 bits resolution
% lower_saturation   = -2^15;   
wn                 = 1000*2*pi;   % MEMS dynamics [rad/s]
% w                  = 12.5*2*pi;   % LPF cut off frequency [rad/s]
zeta               = sqrt(2)/2;   % Damping estimate
% z                  = zeta;

t=1:length(true_w);
for i = 1:length(true_w)
    noise(i,1) = rate_noise_density^2*randn();
end
H = tf([wn^2],[1 2*zeta*wn wn^2]);
sensed_angular_rate = lsim(H,true_w,t) + gyro_bias + noise;

G_rate = sensed_angular_rate;
% for i = 1:length(true_w)
%     G_rate(:,i) = true_w(:,i) + .01*(pi()/180)*randn();
% end

end

