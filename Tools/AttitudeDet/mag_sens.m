function [ b_sens_v ] = mag_sens( b_v, n_h )
%MAG_SENS Magnetic Field Strength Model
%n_hat 3xN unit vector normal to sensor
    n_std_dev = 0; %125e-9; %Noise floor (ADIS16404)
    %n_std_dev = 410e-9; %STM
    
    for i = 1:length(b_v(1,:))
        n_b = n_std_dev*randn();
        b_sens = dot(b_v(:,i),n_h(:,i)) + n_b;
	b_sens_v(1,i) = b_sens;
    end
end

