function [ s_e, sens_num ] = est_sun_sens( I_sun_sens, norm_sun_sens )
% EST_SUN_SENS 
% input: I_sun_sens 1xN
% input: norm_sun_sens Kx3xN
% output s_e 3xN
    I_0 = max(max(abs(I_sun_sens)));%3;%Constant based on I_0 current from Sun Sensor

    for j = 1:length(I_sun_sens(1,:))
	[b, A] = cal_sun_sens(I_sun_sens(:,j), norm_sun_sens);
        A_size = size(A);
        sens_num(:,j) = A_size(1,1);

        if sens_num(:,j) > 2
            s_e(:,j) = (1/I_0)*inv(A'*A)*A'*b;
        else
            s_e(:,j) = [0;0;0];
        end
    end
end

