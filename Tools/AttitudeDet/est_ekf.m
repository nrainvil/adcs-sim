function [ DCM_k_p, P_k_p ] = est_ekf( DCM_k_m, drift_bias, P_km_p, w_in, b_k_in, sun_sens, norm_sun_sens, b_ref_v_in, s_ref_v_in, dt)
%EST_MEKF 
% DCM_k_m - DCM of prior state
% w   - body rates 3x1

	%Sun Sensor
        [S_sens_est_bf, S_sens_num_bf] = est_sun_sens(sun_sens, norm_sun_sens);
        
	%Covariance 
	R_b = [((.44e-6)^2)*eye(3,3)]; %Magnetometer - STM
        R_s = [(.02^2)*eye(3,3)]; %Sun Sensor - Osram 
	R = [R_b,zeros(3,3);zeros(3,3),R_s];
	Q = (2*.025)*eye(4,4); %Process Noise

	%Format input data
	w = -w_in; 
	b_ref_v = b_ref_v_in/norm(b_ref_v_in);
	s_ref_v = s_ref_v_in/norm(s_ref_v_in);
	b_k = b_k_in/norm(b_k_in);
	s_k = S_sens_est_bf/norm(S_sens_est_bf);
	z_ref_v = [b_ref_v;s_ref_v];
	z_k = [b_k;s_k];

	%Propagation
	x_km_p = dcm2q((DCM_k_m));	
	Omega = [      0, w(3,1),-w(2,1), w(1,1); ... 
                 -w(3,1),      0, w(1,1), w(2,1); ...
                  w(2,1),-w(1,1),      0, w(3,1); ...
                 -w(1,1),-w(2,1),-w(3,1),      0];
	F = (1/2)*Omega;
	G = F; %?
	dx = F*x_km_p;
	dP = F*P_km_p+P_km_p*F' + G*Q*G'; 
	x_k_m = dx*dt + x_km_p;
	P_k_m = dP*dt + P_km_p;

%	%Gain - Full
%	H = 2*[x_k_m(1,:)*z_ref_v(1,:)+x_k_m(2,:)*z_ref_v(2,:)+x_k_m(3,:)*z_ref_v(3,:),-x_k_m(2,:)*z_ref_v(1,:)+x_k_m(1,:)*z_ref_v(2,:)+x_k_m(4,:)*z_ref_v(3,:), -x_k_m(3,:)*z_ref_v(1,:)-x_k_m(4,:)*z_ref_v(2,:)+x_k_m(1,:)*z_ref_v(3,:),x_k_m(4,:)*z_ref_v(1,:)-x_k_m(3,:)*z_ref_v(2,:)+x_k_m(2,:)*z_ref_v(3,:);...
%               x_k_m(2,:)*z_ref_v(1,:)-x_k_m(1,:)*z_ref_v(2,:)-x_k_m(4,:)*z_ref_v(3,:),x_k_m(1,:)*z_ref_v(1,:)+x_k_m(2,:)*z_ref_v(2,:)+x_k_m(3,:)*z_ref_v(3,:),x_k_m(4,:)*z_ref_v(1,:)-x_k_m(3,:)*z_ref_v(2,:)+x_k_m(2,:)*z_ref_v(3,:),x_k_m(3,:)*z_ref_v(1,:)+x_k_m(4,:)*z_ref_v(2,:)-x_k_m(1,:)*z_ref_v(3,:); ...
%               x_k_m(3,:)*z_ref_v(1,:)+x_k_m(4,:)*z_ref_v(2,:)-x_k_m(1,:)*z_ref_v(3,:),-x_k_m(4,:)*z_ref_v(1,:)+x_k_m(3,:)*z_ref_v(2,:)-x_k_m(2,:)*z_ref_v(3,:),x_k_m(1,:)*z_ref_v(1,:)+x_k_m(2,:)*z_ref_v(2,:)+x_k_m(3,:)*z_ref_v(3,:),-x_k_m(2,:)*z_ref_v(1,:)+x_k_m(1,:)*z_ref_v(2,:)+x_k_m(4,:)*z_ref_v(3,:); ...
%               x_k_m(1,:)*z_ref_v(4,:)+x_k_m(2,:)*z_ref_v(5,:)+x_k_m(3,:)*z_ref_v(6,:),-x_k_m(2,:)*z_ref_v(4,:)+x_k_m(1,:)*z_ref_v(5,:)+x_k_m(4,:)*z_ref_v(6,:), -x_k_m(3,:)*z_ref_v(4,:)-x_k_m(4,:)*z_ref_v(5,:)+x_k_m(1,:)*z_ref_v(6,:),x_k_m(4,:)*z_ref_v(4,:)-x_k_m(3,:)*z_ref_v(5,:)+x_k_m(2,:)*z_ref_v(6,:);...
%               x_k_m(2,:)*z_ref_v(4,:)-x_k_m(1,:)*z_ref_v(5,:)-x_k_m(4,:)*z_ref_v(6,:),x_k_m(1,:)*z_ref_v(4,:)+x_k_m(2,:)*z_ref_v(5,:)+x_k_m(3,:)*z_ref_v(6,:),x_k_m(4,:)*z_ref_v(4,:)-x_k_m(3,:)*z_ref_v(5,:)+x_k_m(2,:)*z_ref_v(6,:),x_k_m(3,:)*z_ref_v(4,:)+x_k_m(4,:)*z_ref_v(5,:)-x_k_m(1,:)*z_ref_v(6,:); ...
%               x_k_m(3,:)*z_ref_v(4,:)+x_k_m(4,:)*z_ref_v(5,:)-x_k_m(1,:)*z_ref_v(6,:),-x_k_m(4,:)*z_ref_v(4,:)+x_k_m(3,:)*z_ref_v(5,:)-x_k_m(2,:)*z_ref_v(6,:),x_k_m(1,:)*z_ref_v(4,:)+x_k_m(2,:)*z_ref_v(5,:)+x_k_m(3,:)*z_ref_v(6,:),-x_k_m(2,:)*z_ref_v(4,:)+x_k_m(1,:)*z_ref_v(5,:)+x_k_m(4,:)*z_ref_v(6,:)];
%
%	K_k = P_k_m*H'*inv(H*P_k_m*H'+R);
%
%	%Update - Full
%	h_quad = [(x_k_m(1,:)^2-x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(2*x_k_m(1,:)*x_k_m(2,:)-2*x_k_m(3,:)*x_k_m(4,:)),(2*x_k_m(2,:)*x_k_m(4,:)+2*x_k_m(1,:)*x_k_m(3,:)); ...
%		  (2*x_k_m(1,:)*x_k_m(2,:)+2*x_k_m(3,:)*x_k_m(4,:)),(-x_k_m(1,:)^2+x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(-2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)); ...
%		  (2*x_k_m(1,:)*x_k_m(3,:)-2*x_k_m(2,:)*x_k_m(4,:)),(2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)),(-x_k_m(1,:)^2-x_k_m(2,:)^2+x_k_m(3,:)^2+x_k_m(4,:)^2)]; %...
%	h_quad = [h_quad,zeros(3,3); ...
%                  zeros(3,3),h_quad];
%
%        x_quad_k = [h_quad*z_ref_v];
%	x_k_p = x_k_m + K_k*(z_k-x_quad_k);
%	P_k_p = (eye(4,4)-K_k*H)*P_k_m; 

	%Gain - B 
	H_b = 2*[x_k_m(1,:)*b_ref_v(1,:)+x_k_m(2,:)*b_ref_v(2,:)+x_k_m(3,:)*b_ref_v(3,:),-x_k_m(2,:)*b_ref_v(1,:)+x_k_m(1,:)*b_ref_v(2,:)+x_k_m(4,:)*b_ref_v(3,:), -x_k_m(3,:)*b_ref_v(1,:)-x_k_m(4,:)*b_ref_v(2,:)+x_k_m(1,:)*b_ref_v(3,:),x_k_m(4,:)*b_ref_v(1,:)-x_k_m(3,:)*b_ref_v(2,:)+x_k_m(2,:)*b_ref_v(3,:);...
               x_k_m(2,:)*b_ref_v(1,:)-x_k_m(1,:)*b_ref_v(2,:)-x_k_m(4,:)*b_ref_v(3,:),x_k_m(1,:)*b_ref_v(1,:)+x_k_m(2,:)*b_ref_v(2,:)+x_k_m(3,:)*b_ref_v(3,:),x_k_m(4,:)*b_ref_v(1,:)-x_k_m(3,:)*b_ref_v(2,:)+x_k_m(2,:)*b_ref_v(3,:),x_k_m(3,:)*b_ref_v(1,:)+x_k_m(4,:)*b_ref_v(2,:)-x_k_m(1,:)*b_ref_v(3,:); ...
               x_k_m(3,:)*b_ref_v(1,:)+x_k_m(4,:)*b_ref_v(2,:)-x_k_m(1,:)*b_ref_v(3,:),-x_k_m(4,:)*b_ref_v(1,:)+x_k_m(3,:)*b_ref_v(2,:)-x_k_m(2,:)*b_ref_v(3,:),x_k_m(1,:)*b_ref_v(1,:)+x_k_m(2,:)*b_ref_v(2,:)+x_k_m(3,:)*b_ref_v(3,:),-x_k_m(2,:)*b_ref_v(1,:)+x_k_m(1,:)*b_ref_v(2,:)+x_k_m(4,:)*b_ref_v(3,:)];

	K_b_k = (P_k_m*H_b'*inv(H_b*P_k_m*H_b'+R_b)); 

	%Update - B
	h_quad = [(x_k_m(1,:)^2-x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(2*x_k_m(1,:)*x_k_m(2,:)-2*x_k_m(3,:)*x_k_m(4,:)),(2*x_k_m(2,:)*x_k_m(4,:)+2*x_k_m(1,:)*x_k_m(3,:)); ...
		  (2*x_k_m(1,:)*x_k_m(2,:)+2*x_k_m(3,:)*x_k_m(4,:)),(-x_k_m(1,:)^2+x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(-2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)); ...
		  (2*x_k_m(1,:)*x_k_m(3,:)-2*x_k_m(2,:)*x_k_m(4,:)),(2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)),(-x_k_m(1,:)^2-x_k_m(2,:)^2+x_k_m(3,:)^2+x_k_m(4,:)^2)];
        x_quad_k = [h_quad*b_ref_v(1:3,:)];
	x_b_k_p = x_k_m + K_b_k*(b_k-x_quad_k);
	P_b_k_p = (eye(4,4)-K_b_k*H_b)*P_k_m; 
	
	%Gain - S
	H_s = 2*[x_b_k_p(1,:)*s_ref_v(1,:)+x_b_k_p(2,:)*s_ref_v(2,:)+x_b_k_p(3,:)*s_ref_v(3,:),-x_b_k_p(2,:)*s_ref_v(1,:)+x_b_k_p(1,:)*s_ref_v(2,:)+x_b_k_p(4,:)*s_ref_v(3,:), -x_b_k_p(3,:)*s_ref_v(1,:)-x_b_k_p(4,:)*s_ref_v(2,:)+x_b_k_p(1,:)*s_ref_v(3,:),x_b_k_p(4,:)*s_ref_v(1,:)-x_b_k_p(3,:)*s_ref_v(2,:)+x_b_k_p(2,:)*s_ref_v(3,:);...
               x_b_k_p(2,:)*s_ref_v(1,:)-x_b_k_p(1,:)*s_ref_v(2,:)-x_b_k_p(4,:)*s_ref_v(3,:),x_b_k_p(1,:)*s_ref_v(1,:)+x_b_k_p(2,:)*s_ref_v(2,:)+x_b_k_p(3,:)*s_ref_v(3,:),x_b_k_p(4,:)*s_ref_v(1,:)-x_b_k_p(3,:)*s_ref_v(2,:)+x_b_k_p(2,:)*s_ref_v(3,:),x_b_k_p(3,:)*s_ref_v(1,:)+x_b_k_p(4,:)*s_ref_v(2,:)-x_b_k_p(1,:)*s_ref_v(3,:); ...
               x_b_k_p(3,:)*s_ref_v(1,:)+x_b_k_p(4,:)*s_ref_v(2,:)-x_b_k_p(1,:)*s_ref_v(3,:),-x_b_k_p(4,:)*s_ref_v(1,:)+x_b_k_p(3,:)*s_ref_v(2,:)-x_b_k_p(2,:)*s_ref_v(3,:),x_b_k_p(1,:)*s_ref_v(1,:)+x_b_k_p(2,:)*s_ref_v(2,:)+x_b_k_p(3,:)*s_ref_v(3,:),-x_b_k_p(2,:)*s_ref_v(1,:)+x_b_k_p(1,:)*s_ref_v(2,:)+x_b_k_p(4,:)*s_ref_v(3,:)];

	K_s_k = (P_b_k_p*H_s'*inv(H_s*P_b_k_p*H_s'+R_s)); 

	%Update - S
	h_quad = [(x_b_k_p(1,:)^2-x_b_k_p(2,:)^2-x_b_k_p(3,:)^2+x_b_k_p(4,:)^2),(2*x_b_k_p(1,:)*x_b_k_p(2,:)-2*x_b_k_p(3,:)*x_b_k_p(4,:)),(2*x_b_k_p(2,:)*x_b_k_p(4,:)+2*x_b_k_p(1,:)*x_b_k_p(3,:)); ...
		  (2*x_b_k_p(1,:)*x_b_k_p(2,:)+2*x_b_k_p(3,:)*x_b_k_p(4,:)),(-x_b_k_p(1,:)^2+x_b_k_p(2,:)^2-x_b_k_p(3,:)^2+x_b_k_p(4,:)^2),(-2*x_b_k_p(1,:)*x_b_k_p(4,:)+2*x_b_k_p(2,:)*x_b_k_p(3,:)); ...
		  (2*x_b_k_p(1,:)*x_b_k_p(3,:)-2*x_b_k_p(2,:)*x_b_k_p(4,:)),(2*x_b_k_p(1,:)*x_b_k_p(4,:)+2*x_b_k_p(2,:)*x_b_k_p(3,:)),(-x_b_k_p(1,:)^2-x_b_k_p(2,:)^2+x_b_k_p(3,:)^2+x_b_k_p(4,:)^2)];
        x_quad_k = [h_quad*s_ref_v(1:3,:)];
	x_k_p = x_b_k_p + K_s_k*(s_k-x_quad_k);
	P_k_p = (eye(4,4)-K_s_k*H_s)*P_b_k_p; 
	
	%Norm
	x_k_p = x_k_p/norm(x_k_p);
	DCM_k_p = q2dcm(x_k_p);
	
end

