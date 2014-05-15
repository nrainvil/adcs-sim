function [ I, alpha ] = sun_sens( s_hat, n_hat_3d)
%SUN_SENS_COARSE Current model for a coarse sun sensor
%s_hat 3xN unit vectors
%n_hat Kx3xN
    I_0 = 3; % V RAX2 - OSRAMSFH2430
    N = length(s_hat(1,:));
    K = length(n_hat_3d(:,1,1));
    rng('default'); %Reset randn
    n_std_dev = .02; %.005; %.015/3 SolarMEMs A60 Estimate (.02 Osram)  V Noise std-dev
    half_angle = 70; % Deg (FOV) (70 Osram)
    
    for k = 1:K
        n_hat = reshape(n_hat_3d(k,:,:),3,N);
        for n = 1:N
            n_g = n_std_dev*randn(); % 0 mean
            sens_cos = dot(n_hat(:,n),s_hat(:,n));
            alpha(1,n) = acosd(sens_cos);
            if alpha(1,n) > -half_angle  && alpha(1,n) < half_angle 
                I(k,n) = I_0*sens_cos + n_g;
            else
                I(k,n) = 0;
            end 
        end
    end
end

