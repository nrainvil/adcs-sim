function Q_out = Math_Q_Prop(Q_in)
%#eml

Q_out = zeros(size(Q_in));

if (Q_in(4) >= 0)
	Q_out(1) = Q_in(1);
	Q_out(2) = Q_in(2);
	Q_out(3) = Q_in(3);
	Q_out(4) = Q_in(4);
else
	Q_out(1) = -1.0 * Q_in(1);
	Q_out(2) = -1.0 * Q_in(2);
	Q_out(3) = -1.0 * Q_in(3);
	Q_out(4) = -1.0 * Q_in(4);
end
