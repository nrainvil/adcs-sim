function [ alpha_m ] = cart2ss( alpha_v )
%CART2SS 
	alpha_m = [            0, -alpha_v(3,1),  alpha_v(2,1); ...
                    alpha_v(3,1),             0, -alpha_v(1,1); ...
                   -alpha_v(2,1),  alpha_v(1,1),             0];

end

