function [B_mag_sens_rs, I_sun_sens_rs, G_rate_sens_rs, B_ECI_rs, Sun_ECI_rs, Time_rs] = adc_sim(fs_sens, B_mag_sens, I_sun_sens, G_rate_sens, Orbit)
%ADC_SIM 
	%Convert to Binary
	%Magnetometer 16-bit 2's Complment
	%mag_sens = .1e-6; %Freescale
	mag_bit_width = 16;
	%mag_bit_width = 8;
	max_mag = max(max(abs(B_mag_sens)));
	mag_sens = max_mag/(2^(mag_bit_width-1)); %FIX
	B_mag_sens_disc = twos_dec2bin(B_mag_sens./mag_sens,mag_bit_width);

	%Sun Sensors 12-bit 2's Complement
	sun_bit_width = 12;
	%sun_bit_width = 8;
	max_sun = max(max(abs(I_sun_sens)));
	sun_sens = max_sun/(2^(sun_bit_width-1)); %FIX
	I_sun_sens_disc = twos_dec2bin(I_sun_sens./sun_sens,sun_bit_width);

	%Convert to Decimal
	B_mag_sens_conv = twos_bin2dec(B_mag_sens_disc,mag_bit_width);
	I_sun_sens_conv = twos_bin2dec(I_sun_sens_disc,sun_bit_width);
	%I_sun_sens_conv = I_sun_sens;

	%Shared Sample Rate
	Sun_ECI_ts = timeseries(Orbit.Sun_ECI,Orbit.Time,'name','I_eci');
	I_sun_sens_ts = timeseries(I_sun_sens_conv,Orbit.Time,'name','I_bf');
	B_ECI_ts = timeseries(Orbit.B_ECI,Orbit.Time,'name','B_eci');
	B_mag_sens_ts = timeseries(B_mag_sens_conv,Orbit.Time,'name','B_bf');
	G_rate_sens_ts = timeseries(G_rate_sens,Orbit.Time,'name','G_bf');

	Time_rs = Orbit.Time(1):(1/fs_sens):(Orbit.Time(length(Orbit.Time)));
	Sun_ECI_rs = resample(Sun_ECI_ts,Time_rs);
	I_sun_sens_rs = resample(I_sun_sens_ts,Time_rs);
	B_ECI_rs = resample(B_ECI_ts,Time_rs);
	B_mag_sens_rs = resample(B_mag_sens_ts,Time_rs);
	G_rate_sens_rs = resample(G_rate_sens_ts,Time_rs);

	

end

