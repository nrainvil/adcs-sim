%Load mag_test.txt 10Hz Freescale Mag3110FCR1CT-ND
%Expected .25 uT Noise RMS
b_v(1,:) = xaxis(2:length(xaxis))';
b_v(2,:) = yaxis(2:length(yaxis))';
b_v(3,:) = zaxis(2:length(zaxis))';

t = 0:.1:((length(b_x)-1)/10);

b_v_norm(1,:) = b_v(1,:) - mean(b_v(1,:));
b_v_norm(2,:) = b_v(2,:) - mean(b_v(2,:));
b_v_norm(3,:) = b_v(3,:) - mean(b_v(3,:));

b_v_rms_err(1,:) = sqrt((1/length(b_v_norm(1,:)))*sum((b_v_norm(1,:).^2)));
b_v_rms_err(2,:) = sqrt((1/length(b_v_norm(2,:)))*sum((b_v_norm(2,:).^2)));
b_v_rms_err(3,:) = sqrt((1/length(b_v_norm(3,:)))*sum((b_v_norm(3,:).^2)));

fprintf('X Axis RMS Error = %.3f uT\n',b_v_rms_err(1,:));
fprintf('Y Axis RMS Error = %.3f uT\n',b_v_rms_err(2,:));
fprintf('Z Axis RMS Error = %.3f uT\n',b_v_rms_err(3,:));

figure;
hold on;
plot(t,b_v(1,:),'r');
plot(t,b_v(2,:),'g');
plot(t,b_v(3,:),'b');
hold off;
xlabel('Time (s)');
ylabel('Field Strength (uT)');

