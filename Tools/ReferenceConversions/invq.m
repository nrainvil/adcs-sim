function [ q_out ] = invq( q )
%INVQ Summary of this function goes here
%   Detailed explanation goes here
	q_out(1,:) = -q(1,:);
	q_out(2,:) = -q(2,:);
	q_out(3,:) = -q(3,:);
	q_out(4,:) = q(4,:);


end

