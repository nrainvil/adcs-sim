function [ b_sens_v ] = mag_sens_3axis( b_v, x_h, z_h)
%MAG_SENS_3AXIS 3-Axis Magnetometer Model
%b_v - magnetic field vector, x_h vector in X-Axis of sensor.
    n_std_dev = 125e-9; %Noise floor (ADIS16404)
    %n_std_dev = 410e-9; %STM

    y_h = cross(z_h,x_h);

    scale = [1;1;1]; %scale factor
    offset = [0;0;0]; %offset
    nonorth = [0;0;0]; %non-orthogonality angle (degrees)

    Bx_v = sum((b_v.*x_h));
    By_v = sum((b_v.*y_h));
    Bz_v = sum((b_v.*z_h));

    b_sens_v = zeros(3,length(b_v(1,:)));
    b_sens_c = num2cell(b_sens_v,1);

    parfor i = 1:length(b_v(1,:))
        b_sens_c{i}(1) = scale(1,1)*Bx_v(:,i) + offset(1,1) + n_std_dev*randn();
	b_sens_c{i}(2) = scale(2,1)*(By_v(:,i)*cosd(nonorth(1,1)+Bx_v(:,i)*sind(nonorth(1,1)))) + offset(2,1) + n_std_dev*randn();
	b_sens_c{i}(3) = scale(3,1)*(Bx_v(:,i)*sind(nonorth(2,1))+By_v(:,i)*sind(nonorth(3,1))*cosd(nonorth(2,1))+Bz_v(:,i)*cosd(nonorth(3,1))*cosd(nonorth(2,1))) + offset(3,1) + n_std_dev*randn();
    end
    b_sens_v = [b_sens_c{:}];
end

