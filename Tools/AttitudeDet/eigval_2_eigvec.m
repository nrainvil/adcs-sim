function [ eig_vec ] = eigval_2_eigvec(eig_val,S,Z,sigma)
%EIGVAL_2_EIGVEC 

p = inv((eig_val+sigma)*eye(3,3)-S)*Z;
eig_vec = (1/sqrt(1+p'*p))*vertcat(p,1);

end

