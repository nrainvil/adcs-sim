function [ DCM_k_p, P_k_p ] = est_ekf( DCM_k_m, drift_bias, P_km_p, w_in, z_k_in, ref_v_in, dt)
%EST_MEKF 
% DCM_k_m - DCM of prior state
% w   - body rates 3x1
	%Covariance 
	R_b = [((.44e-6)^2)*eye(3,3)]; %Magnetometer - STM
        R_s = .02^2; %Sun Sensor - Osram 
	Q = (2*.025)*eye(4,4); %Process Noise

	%Format input data
	w = -w_in; 
	ref_v = ref_v_in/norm(ref_v_in);
	z_k = z_k_in/norm(z_k_in);

	%Propagation
	x_km_p = dcm2q((DCM_k_m));	
	Omega = [      0, w(3,1),-w(2,1), w(1,1); ... 
                 -w(3,1),      0, w(1,1), w(2,1); ...
                  w(2,1),-w(1,1),      0, w(3,1); ...
                 -w(1,1),-w(2,1),-w(3,1),      0];
	F = (1/2)*Omega;
	G = F; %?
	dx = F*x_km_p;
	dP = F*P_km_p+P_km_p*F' + G*Q+G'; 
	x_k_m = dx*dt + x_km_p;
	P_k_m = dP*dt + P_km_p;

	%Gain - B 
	H_b = 2*[x_k_m(1,:)*ref_v(1,:)+x_k_m(2,:)*ref_v(2,:)+x_k_m(3,:)*ref_v(3,:),-x_k_m(2,:)*ref_v(1,:)+x_k_m(1,:)*ref_v(2,:)+x_k_m(4,:)*ref_v(3,:), -x_k_m(3,:)*ref_v(1,:)-x_k_m(4,:)*ref_v(2,:)+x_k_m(1,:)*ref_v(3,:),x_k_m(4,:)*ref_v(1,:)-x_k_m(3,:)*ref_v(2,:)+x_k_m(2,:)*ref_v(3,:);...
               x_k_m(2,:)*ref_v(1,:)-x_k_m(1,:)*ref_v(2,:)-x_k_m(4,:)*ref_v(3,:),x_k_m(1,:)*ref_v(1,:)+x_k_m(2,:)*ref_v(2,:)+x_k_m(3,:)*ref_v(3,:),x_k_m(4,:)*ref_v(1,:)-x_k_m(3,:)*ref_v(2,:)+x_k_m(2,:)*ref_v(3,:),x_k_m(3,:)*ref_v(1,:)+x_k_m(4,:)*ref_v(2,:)-x_k_m(1,:)*ref_v(3,:); ...
               x_k_m(3,:)*ref_v(1,:)+x_k_m(4,:)*ref_v(2,:)-x_k_m(1,:)*ref_v(3,:),-x_k_m(4,:)*ref_v(1,:)+x_k_m(3,:)*ref_v(2,:)-x_k_m(2,:)*ref_v(3,:),x_k_m(1,:)*ref_v(1,:)+x_k_m(2,:)*ref_v(2,:)+x_k_m(3,:)*ref_v(3,:),-x_k_m(2,:)*ref_v(1,:)+x_k_m(1,:)*ref_v(2,:)+x_k_m(4,:)*ref_v(3,:)];

	K_b_k = P_k_m*H_b'*inv(H_b*P_k_m*H_b'+R); 

	%Update - B
	h_quad = [(x_k_m(1,:)^2-x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(2*x_k_m(1,:)*x_k_m(2,:)-2*x_k_m(3,:)*x_k_m(4,:)),(2*x_k_m(2,:)*x_k_m(4,:)+2*x_k_m(1,:)*x_k_m(3,:)); ...
		  (2*x_k_m(1,:)*x_k_m(2,:)+2*x_k_m(3,:)*x_k_m(4,:)),(-x_k_m(1,:)^2+x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(-2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)); ...
		  (2*x_k_m(1,:)*x_k_m(3,:)-2*x_k_m(2,:)*x_k_m(4,:)),(2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)),(-x_k_m(1,:)^2-x_k_m(2,:)^2+x_k_m(3,:)^2+x_k_m(4,:)^2)];
        x_quad_k = [h_quad*ref_v(1:3,:)];
	x_b_k_p = x_k_m + K_b_k*(z_k-x_quad_k);
	P_b_k_p = (eye(4,4)-K_b_k*H_b)*P_k_m; 

	%Norm
	x_k_p = x_b_k_p; %TEMP
	x_k_p = x_k_p/norm(x_k_p);
	DCM_k_p = q2dcm(x_k_p);
	
	P_k_p = P_b_k_p; %TEMP

end

