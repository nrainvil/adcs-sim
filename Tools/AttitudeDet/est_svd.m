function [ R ] = est_svd(b_k, eci_k)
%EST_SVD 
% input: b_k = 3xN input vector, body frame
% input: eci_k = 3xN input vector, eci frame
% output: R = rotation matrix, body to eci

B = zeros(3,3);
for k=1:length(b_k(1,:))
    B = B + eci_k(:,k)*b_k(:,k)';
end

[U,S,V] = svd(B);
M = diag([1,1,det(U)*det(V)]);
R = U*M*V';

end

