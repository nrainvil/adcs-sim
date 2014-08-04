function [ I_sun_sens_out, norm_sun_sens_out] = cal_sun_sens( I_sun_sens, norm_sun_sens )
%CAL_SUN_SENS
    I_0 = max(max(abs(I_sun_sens)));%3;%Constant based on I_0 current from Sun Sensor
    I_co = I_0/10;%Cuttoff current, all sensors below this value are ignored
    I_sun_sens_out = [];
    norm_sun_sens_out = [];

    for i = 1:length(I_sun_sens(:,1))
	curr_sens = norm_sun_sens(i,:);
	curr_sens_diff = ones(length(norm_sun_sens(:,1)),1)*curr_sens + norm_sun_sens;
	norm_sens_diff = sum(abs(curr_sens_diff)')';
	[~,rank] = sort(norm_sens_diff);
	opposite_sens_loc = rank(1,:);
	%Test opposite sensor
	if (abs(I_sun_sens(i,:)) > I_co ) & (I_sun_sens(i,:) > I_sun_sens(opposite_sens_loc,:)) 
	   I_sun_sens_out = [I_sun_sens_out;I_sun_sens(i,:)];
	   norm_sun_sens_out = [norm_sun_sens_out;norm_sun_sens(i,:)];
	end
    end
end

