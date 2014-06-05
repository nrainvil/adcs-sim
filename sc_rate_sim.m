function [ G_rate_sens] = sc_sens_sim(R_dot_eci_body)

%3x Rate Gyros
G_rate_sens(1,:) = rate_sens(reshape(R_dot_eci_body(3,2,:),1,length(R_dot_eci_body(3,2,:))));
G_rate_sens(2,:) = rate_sens(reshape(R_dot_eci_body(1,3,:),1,length(R_dot_eci_body(1,3,:))));
G_rate_sens(3,:) = rate_sens(reshape(R_dot_eci_body(2,1,:),1,length(R_dot_eci_body(2,1,:))));

end
