function [ dec_v ] = twos_bin2dec( bin_v, N)
%TWOS_BIN2DEC 
if length(size(bin_v)) > 2
	for i=1:length(bin_v(1,:,1))
		for j=1:length(bin_v(:,1,1))
			bin = reshape(bin_v(j,i,:),1,N);
			if strcmp(bin(1,1),'1')
				dec = bin2dec(bin)-2^length(bin);	
			else
				dec = bin2dec(bin);
			end
			dec_v(j,i) = dec;
		end
		%dec_v(:,i) = dec_v_j;
	end
else
	bin = bin_v;
	if strcmp(bin(1,1),'1')
		dec = bin2dec(bin)-2^length(bin);	
	else
		dec = bin2dec(bin);
	end
	dec_v = dec;
end

end

