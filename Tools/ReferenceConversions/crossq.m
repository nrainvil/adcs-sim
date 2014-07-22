function [ q_out ] = crossq( q, p)
%CROSSQ Multiply two quaternions
	q_v = q(1:3,:);
	p_v = p(1:3,:);
	q_w = q(4,:);
	p_w = p(4,:);

	q_out_v = cross(q_v,p_v) + q_w.*p_v + p_w.*q_v;
	q_out_w = q_w*p_w-dot(q_v,p_v);
	q_out = vertcat(q_out_v,q_out_w); 
end
