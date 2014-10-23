function [ R ] = est_quest( b_k, eci_k )
%EST_QUEST
% input: b_k = 3xN input vector, body frame
% input: eci_k = 3xN input vector, eci frame
% output: R = rotation matrix, body to eci

B = zeros(3,3);
a_i = (1/length(b_k(1,:))); %Normalized Weight
for k=1:length(b_k(1,:))
    w_v = b_k(:,k)/norm(b_k(:,k));
    v_v = eci_k(:,k)/norm(eci_k(:,k));
    B = B + a_i*w_v*v_v';
end

sigma = trace(B);
S = B + B';
Z = [B(2,3)-B(3,2); ...
     B(3,1)-B(1,3); ...
     B(1,2)-B(2,1)];

%%Assume eigenvalue is 1
eig_val = 1; 
q = eigval_2_eigvec(eig_val,S,Z,sigma);

%1st iteration Newton
%d_eig = .01*eig_val;
%q_pde = eigval_2_eigvec(eig_val+d_eig/2,S,Z,sigma);
%q_mde = eigval_2_eigvec(eig_val-d_eig/2,S,Z,sigma);
%eig_val_p = eig_val - sum(K*q-eig_val*q)./((sum(K*q_pde-(eig_val+d_eig/2)*q_pde)-sum(K*q_mde-(eig_val-d_eig/2)*q_mde))/d_eig);
%eig_val = eig_val_p;
%q = eigval_2_eigvec(eig_val,S,Z,sigma);

R = q2dcm(q)'; %Transpose!!!

end

