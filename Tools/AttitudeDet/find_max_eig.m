function [ eig_max ] = find_max_eig( v, w )
%FIND_MAX_EIG 
%Even weight, 2 column vectors

a_i = .5;
vec_cos = dot(v(:,1),v(:,2))*dot(w(:,1),w(:,2))+norm(cross(v(:,1),v(:,2)))*norm(cross(w(:,1),w(:,2)));
eig_max = sqrt(a_i^2 + 2*a_i*a_i*vec_cos+a_i^2);

%k = trace(S');
%del = det(S);
%sig = (1/2)*trace(S);
%a = sig^2 - k;
%b = sig^2 + Z'*Z;
%c = del+Z'*S*Z;
%d = Z'*(S'*S)*Z;



end

