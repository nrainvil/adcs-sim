function [ I, alpha ] = sun_sens( s_hat, n_hat)
%SUN_SENS_COARSE Current model for a coarse sun sensor
%s_hat, n_hat 3xN unit vectors
    I_0 = 3; % V RAX2 - OSRAMSFH2430
    rng('default'); %Reset randn
    n_std_dev = .02; % V Noise std-dev
    
    for i = 1:length(s_hat(1,:))
        n_g = n_std_dev*randn(); % 0 mean
        sens_cos = dot(n_hat(:,i),s_hat(:,i));
        alpha(1,i) = acosd(sens_cos);
        if alpha(1,i) > -70  && alpha(1,i) < 70
            I(:,i) = I_0*sens_cos + n_g;
        else
            I(:,i) = 0;
        end
    end
end

