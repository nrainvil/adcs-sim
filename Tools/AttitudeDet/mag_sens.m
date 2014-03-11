function [ b_sens_v ] = mag_sens( b_v, n_h )
%MAG_SENS Magnetic Field Strength Model
%n_hat 3xN unit vector normal to sensor
    n_std_dev = 111e-9; %Noise floor (ADIS16404)
    
    for i = 1:length(b_v)
        n_b = n_std_dev*randn();
        b_sens_v(1,i) = dot(b_v(:,i),n_h(:,i)) + n_b;
    end
end
