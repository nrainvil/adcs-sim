function [ R ] = est_quest( b_k, eci_k )
%EST_QUEST
% input: b_k = 3xN input vector, body frame
% input: eci_k = 3xN input vector, eci frame
% output: R = rotation matrix, body to eci

B = zeros(3,3);
for k=1:length(b_k(1,:))
    B = B + (1/length(b_k(1,:)))*b_k(:,k)*eci_k(:,k)'; %Normalized weight
end
S = B + B';
Z = [B(2,3)-B(3,2); ...
     B(3,1)-B(1,3); ...
     B(1,2)-B(2,1)];
sigma = trace(B);

%%TEMP
K = vertcat(horzcat(S-sigma*eye(3,3), Z),horzcat(Z',sigma));

%Find max eigenvalue of K
[V,D] = eig(K);
[m0,l0] = max(D);
[m,l] = max(m0);
eig_val_true = D(l0(l),l);

%%TEMP
eig_val = 1;
d_eig = .01*eig_val;
q = eigval_2_eigvec(eig_val,S,Z,sigma);
%1st iteration Newton
q_pde = eigval_2_eigvec(eig_val+d_eig/2,S,Z,sigma);
q_mde = eigval_2_eigvec(eig_val-d_eig/2,S,Z,sigma);
eig_val_p = eig_val - sum(K*q-eig_val*q)./((sum(K*q_pde-(eig_val+d_eig/2)*q_pde)-sum(K*q_mde-(eig_val-d_eig/2)*q_mde))/d_eig);
eig_val = eig_val_p;
q = eigval_2_eigvec(eig_val,S,Z,sigma);
%2nd iteration Newton
q_pde = eigval_2_eigvec(eig_val+d_eig/2,S,Z,sigma);
q_mde = eigval_2_eigvec(eig_val-d_eig/2,S,Z,sigma);
eig_val_p = eig_val - sum(K*q-eig_val*q)./((sum(K*q_pde-(eig_val+d_eig/2)*q_pde)-sum(K*q_mde-(eig_val-d_eig/2)*q_mde))/d_eig);
eig_val = eig_val_p;
q = eigval_2_eigvec(eig_val,S,Z,sigma);


eig_vals = linspace(1-d_eig/2,1+d_eig/2);
for i = 1:length(eig_vals)
    q_p = eigval_2_eigvec(eig_vals(i),S,Z,sigma);
    eval_v_full = (K*q_p - eig_vals(i)*q_p);
    eval_v(i) = sum(eval_v_full(1:4));
end
figure;plot(eig_vals,eval_v);


R = q2dcm(q)'; %Transpose???

end

