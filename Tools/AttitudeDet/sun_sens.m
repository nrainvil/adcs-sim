function [ I, alpha ] = sun_sens( s_hat, a_v, n_hat_3d)
%SUN_SENS_COARSE Current model for a coarse sun sensor
%s_hat 3xN unit vectors
%n_hat Kx3xN
    CONST.AM0 = 1361;
    I_0 = 3; % V RAX2 - OSRAMSFH2430
    N = length(s_hat(1,:));
    K = length(n_hat_3d(:,1,1));
    %L = length(a_v(1,:));
    rng('default'); %Reset randn
    n_std_dev = 0; %.02; %.005; %.015/3 SolarMEMs A60 Estimate (.02 Osram)  V Noise std-dev
    half_angle = 70; % Deg (FOV) (70 Osram)
    
    for k = 1:K
        n_hat = reshape(n_hat_3d(k,:,:),3,N);
        for n = 1:N
            n_g = n_std_dev*randn(); % 0 mean

	    %Sun current
            sens_cos = dot(n_hat(:,n),s_hat(:,n));
            alpha(1,n) = acosd(sens_cos);
            if alpha(1,n) > -half_angle  && alpha(1,n) < half_angle 
                I_sun = I_0*sens_cos;
            else
                I_sun = 0;
            end 

	    %Albedo current
	    I_alb = 0; %Init
	    albedo_v = (I_0/CONST.AM0)*a_v;
	    albedo_hat_v = albedo_v./(ones(3,1)*sqrt(sum(albedo_v.^2)));
	    albedo_hat_v(isnan(albedo_hat_v)) = 0;
            sens_cos_v = sum((n_hat(:,n)*ones(1,length(albedo_hat_v(1,:)))).*albedo_hat_v);
	    alpha_mask_v = acosd(sens_cos_v);
            alpha_mask_v(alpha_mask_v < -half_angle | alpha_mask_v > half_angle) = 0;
            alpha_mask_v(alpha_mask_v ~= 0) = 1;
	    I_alb = alpha_mask_v.*sum((n_hat(:,n)*ones(1,length(albedo_hat_v(1,:)))).*albedo_v);

	    I(k,n) = I_sun + n_g + sum(I_alb);

%	    for l = 1:L
%		    albedo_hat = (albedo_v(:,l))/norm(albedo_v(:,l));
%		    sens_cos = dot(n_hat(:,n),albedo_hat);
%		    alpha(1,n) = acosd(sens_cos);
%		    if alpha(1,n) > -half_angle  && alpha(1,n) < half_angle 
%			I_alb = dot(n_hat(:,n),albedo_v(:,n)) + I_alb;
%		    end 
%	    end
        end
    end

end

