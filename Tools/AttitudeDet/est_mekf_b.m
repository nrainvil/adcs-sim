function [ DCM_k_p, P_k_p ] = est_mekf_b( DCM_k_m, drift_bias, P_km_p, w_in, z_k_in, ref_v_in, dt)
%EST_MEKF 
% DCM_k_m - DCM of prior state
% w   - body rates 3x1
	%Covariance 
	R = [((.44e-6)^2)*eye(3,3),zeros(3,3); ... %STM
             zeros(3,3), (.02^2)*eye(3,3)]; %Osram 
	Q = (2*.025)*eye(4,4);

	%Format input data
	w = -w_in; 
	if length(ref_v_in(1,:))>1
		ref_v = vertcat(ref_v_in(:,1)/norm(ref_v_in(:,1)),ref_v_in(:,2)/norm(ref_v_in(:,2)));
		z_k = vertcat(z_k_in(:,1)/norm(z_k_in(:,1)),z_k_in(:,2)/norm(z_k_in(:,2)));
	else
		ref_v = ref_v_in/norm(ref_v_in);
		z_k = z_k_in/norm(z_k_in);
	end

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

	%Gain 
	H = 2*[x_k_m(1,:)*ref_v(1,:)+x_k_m(2,:)*ref_v(2,:)+x_k_m(3,:)*ref_v(3,:),-x_k_m(2,:)*ref_v(1,:)+x_k_m(1,:)*ref_v(2,:)+x_k_m(4,:)*ref_v(3,:), -x_k_m(3,:)*ref_v(1,:)-x_k_m(4,:)*ref_v(2,:)+x_k_m(1,:)*ref_v(3,:),x_k_m(4,:)*ref_v(1,:)-x_k_m(3,:)*ref_v(2,:)+x_k_m(2,:)*ref_v(3,:);...
               x_k_m(2,:)*ref_v(1,:)-x_k_m(1,:)*ref_v(2,:)-x_k_m(4,:)*ref_v(3,:),x_k_m(1,:)*ref_v(1,:)+x_k_m(2,:)*ref_v(2,:)+x_k_m(3,:)*ref_v(3,:),x_k_m(4,:)*ref_v(1,:)-x_k_m(3,:)*ref_v(2,:)+x_k_m(2,:)*ref_v(3,:),x_k_m(3,:)*ref_v(1,:)+x_k_m(4,:)*ref_v(2,:)-x_k_m(1,:)*ref_v(3,:); ...
               x_k_m(3,:)*ref_v(1,:)+x_k_m(4,:)*ref_v(2,:)-x_k_m(1,:)*ref_v(3,:),-x_k_m(4,:)*ref_v(1,:)+x_k_m(3,:)*ref_v(2,:)-x_k_m(2,:)*ref_v(3,:),x_k_m(1,:)*ref_v(1,:)+x_k_m(2,:)*ref_v(2,:)+x_k_m(3,:)*ref_v(3,:),-x_k_m(2,:)*ref_v(1,:)+x_k_m(1,:)*ref_v(2,:)+x_k_m(4,:)*ref_v(3,:); ...
	       x_k_m(1,:)*ref_v(4,:)+x_k_m(2,:)*ref_v(5,:)+x_k_m(3,:)*ref_v(6,:),-x_k_m(2,:)*ref_v(4,:)+x_k_m(1,:)*ref_v(5,:)+x_k_m(4,:)*ref_v(6,:), -x_k_m(3,:)*ref_v(4,:)-x_k_m(4,:)*ref_v(5,:)+x_k_m(1,:)*ref_v(6,:),x_k_m(4,:)*ref_v(4,:)-x_k_m(3,:)*ref_v(5,:)+x_k_m(2,:)*ref_v(6,:);...
               x_k_m(2,:)*ref_v(4,:)-x_k_m(1,:)*ref_v(5,:)-x_k_m(4,:)*ref_v(6,:),x_k_m(1,:)*ref_v(4,:)+x_k_m(2,:)*ref_v(5,:)+x_k_m(3,:)*ref_v(6,:),x_k_m(4,:)*ref_v(4,:)-x_k_m(3,:)*ref_v(5,:)+x_k_m(2,:)*ref_v(6,:),x_k_m(3,:)*ref_v(4,:)+x_k_m(4,:)*ref_v(5,:)-x_k_m(1,:)*ref_v(6,:); ...
               x_k_m(3,:)*ref_v(4,:)+x_k_m(4,:)*ref_v(5,:)-x_k_m(1,:)*ref_v(6,:),-x_k_m(4,:)*ref_v(4,:)+x_k_m(3,:)*ref_v(5,:)-x_k_m(2,:)*ref_v(6,:),x_k_m(1,:)*ref_v(4,:)+x_k_m(2,:)*ref_v(5,:)+x_k_m(3,:)*ref_v(6,:),-x_k_m(2,:)*ref_v(4,:)+x_k_m(1,:)*ref_v(5,:)+x_k_m(4,:)*ref_v(6,:);];

	K_k = P_k_m*H'*inv(H*P_k_m*H'+R); 
	K_k = 0;

	%Update
	h_quad = [(x_k_m(1,:)^2-x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(2*x_k_m(1,:)*x_k_m(2,:)-2*x_k_m(3,:)*x_k_m(4,:)),(2*x_k_m(2,:)*x_k_m(4,:)+2*x_k_m(1,:)*x_k_m(3,:)); ...
		  (2*x_k_m(1,:)*x_k_m(2,:)+2*x_k_m(3,:)*x_k_m(4,:)),(-x_k_m(1,:)^2+x_k_m(2,:)^2-x_k_m(3,:)^2+x_k_m(4,:)^2),(-2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)); ...
		  (2*x_k_m(1,:)*x_k_m(3,:)-2*x_k_m(2,:)*x_k_m(4,:)),(2*x_k_m(1,:)*x_k_m(4,:)+2*x_k_m(2,:)*x_k_m(3,:)),(-x_k_m(1,:)^2-x_k_m(2,:)^2+x_k_m(3,:)^2+x_k_m(4,:)^2)];
        x_quad_k = [h_quad*ref_v(1:3,:);h_quad*ref_v(4:6,:)];
	x_k_p = x_k_m + K_k*(z_k-x_quad_k);
	P_k_p = (eye(4,4)-K_k*H)*P_k_m; 

	%Norm
	x_k_p = x_k_p/norm(x_k_p);
	DCM_k_p = q2dcm(x_k_p);

end

