function [ DCM_k_p, P_k_p ] = est_mekf( DCM_k_m, drift_bias, P_km_p, w_in, z_k_in, ref_v, R_k, dt)
%EST_MEKF 
% DCM_k_m - DCM of prior state
% w   - body rates 3x1
%x is error state
	%Initialize Measurement noise covariance
	%R_k = eye(6,6); %Measurement Noise Covariance
	B_k_m = [0;0;0]; %FIX

	%Format input data
	w = -w_in; 
	w_k_p = w + B_k_m; %FIX
	z_k = [];
	for i=1:length(z_k_in(1,:))
		z_k=vertcat(z_k,z_k_in(:,i));
	end

	%State Estimate Extrapolation
	q_km_p = dcm2q(DCM_k_m);
	Omega = [      0, w(3,1),-w(2,1), w(1,1); ... 
                 -w(3,1),      0, w(1,1), w(2,1); ...
                  w(2,1),-w(1,1),      0, w(3,1); ...
                 -w(1,1),-w(2,1),-w(3,1),      0];
	d_q = ((1/2)*Omega*q_km_p)*dt;
	d_alpha_k_m = d_q(1:3); 
	x_k_m = vertcat(d_alpha_k_m,B_k_m); %State variable: delta q, drift bias

	%Measurement sensitivity matrix
	H_k = [];
	for i=1:length(ref_v(1,:))
		H_k_i = horzcat(DCM_k_m*cart2ss(ref_v(:,i)),zeros(3,3));
		H_k = vertcat(H_k,H_k_i);
	end

	%Covariance Estimate Extrapolation
	P_k_m = P_km_p; % Phi_km*P_km_p*Phi_km' + Q_km; %Covariance

	%Filter Gain
	K_k = P_k_m*H_k'*inv(H_k*P_k_m*H_k'+R_k);
	K_k = 0*K_k; %TEMP - Forces body rate only estimation

	%State Estimate Update
	x_k_p = x_k_m + K_k*(z_k-H_k*x_k_m);

	%Covariance Estimate Update
	P_k_p = (eye(6,6)-K_k*H_k)*P_k_m;

	%Attitude update
	psi_k_p = (sin((1/2)*norm(w_k_p)*dt)*w_k_p)/norm(w_k_p);
	Omega_bar_km_p = [cos((1/2)*norm(w_k_p)*dt)*eye(3,3)-cart2ss(psi_k_p), psi_k_p; ...
                          -psi_k_p', cos((1/2)*norm(w_k_p)*dt)];
	q_k_m = Omega_bar_km_p*q_km_p;

	xi_km_p = q_km_p(4,:)*eye(3,3)+[0, -q_km_p(3,:), q_km_p(2,:); ...
                                        q_km_p(3,:), 0, -q_km_p(1,:); ...
                                        -q_km_p(2,:), q_km_p(1,:), 0;];
	q_k_p_un = q_k_m + vertcat((1/2)*xi_km_p*x_k_p(1:3,:),0);
	q_k_p = q_k_p_un/norm(q_k_p_un);
	DCM_k_p = q2dcm(q_k_p);
 
	%Bias update
	B_k_p = B_k_m + x_k_p(4:6);
end

