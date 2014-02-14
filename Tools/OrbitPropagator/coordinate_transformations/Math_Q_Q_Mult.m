function Q_out = Math_Q_Q_Mult(Q_in1, Q_in2)
%#eml

Q_out = zeros(size(Q_in1));

Q_out(1) = Q_in1(4)*Q_in2(1) - Q_in1(3)*Q_in2(2) + Q_in1(2)*Q_in2(3) + Q_in1(1)*Q_in2(4);
Q_out(2) = Q_in1(4)*Q_in2(2) + Q_in1(3)*Q_in2(1) + Q_in1(2)*Q_in2(4) - Q_in1(1)*Q_in2(3);
Q_out(3) = Q_in1(4)*Q_in2(3) + Q_in1(3)*Q_in2(4) - Q_in1(2)*Q_in2(1) + Q_in1(1)*Q_in2(2);
Q_out(4) = Q_in1(4)*Q_in2(4) - Q_in1(3)*Q_in2(3) - Q_in1(2)*Q_in2(2) - Q_in1(1)*Q_in2(1);

Q_out = Math_Q_Prop(Q_out);