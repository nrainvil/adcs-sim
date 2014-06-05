function [ bin_v ] = twos_dec2bin( dec_v,N )
%TWOS_DEC2BIN Summary of this function goes here
%   Detailed explanation goes here
parfor (i=1:length(dec_v(:,1)))
	bin_temp = repmat(char(0),length(dec_v(1,:)),N);
	for j=1:length(dec_v(1,:))
		dec = dec_v(i,j);
		if (dec >= 0)
			bin = dec2bin(dec,N);
			if length(bin) > N
				bin = strcat(char(48),repmat(char(49),1,N-1));
			end
		else
			if (2^N/2) <= -dec
				bin = dec2bin(2^N/2,N)
			else
				bin = dec2bin(2^N+dec,N);
			end
		end
		bin_temp(j,:) = bin;
	end
	bin_v(i,:,:) = bin_temp;
end

end

