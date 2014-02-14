function [phi e_hat] = q2euler(q)
% input: beta = [beta0 beta1 beta2 beta3]

phi = acos(q(4,:))*2;
e_hat(1:3,:) = q(1:3,:)./sin(phi(1,:)/2);
% e_hat(1,:) = q(1,:)./sin(phi(1,:)/2);
% e_hat(2,:) = q(2,:)./sin(phi(1,:)/2);
% e_hat(3,:) = q(3,:)./sin(phi(1,:)/2);

