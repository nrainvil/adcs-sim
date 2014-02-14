function Q_out = Math_Qinv_Q_Mult(Q_in1, Q_in2)
%#eml
Q_temp    = zeros(size(Q_in1));
Q_temp(1) = -Q_in1(1);
Q_temp(2) = -Q_in1(2);
Q_temp(3) = -Q_in1(3);
Q_temp(4) =  Q_in1(4);

Q_out = Math_Q_Q_Mult(Q_temp, Q_in2);