function [ R ] = est_q( b_k, eci_k)
%EST_Q 
% input: b_k = 3xN input vector, body frame
% input: eci_k = 3xN input vector, eci frame
% output: R = rotation matrix, body to eci

B = zeros(3,3);
for k=1:length(b_k(1,:))
    B = B + (1/length(b_k(1,:)))*b_k(:,k)*eci_k(:,k)';
end
S = B + B';
Z = [B(2,3)-B(3,2); ...
     B(3,1)-B(1,3); ...
     B(1,2)-B(2,1)];
sigma = trace(B);

K = vertcat(horzcat(S-sigma*eye(3,3), Z),horzcat(Z',sigma));

%Find max eigenvalue of K
[V,D] = eig(K);
[m,l] = max(max(D));

%Eigenvector is the quaternion
R = q2dcm(V(:,l))'; %Transpose?

end

