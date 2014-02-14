function V_out = Math_Q_Vec_Qinv_Mult(V_in, Q_in)
%#eml

Q_temp = [0;0;0;0];

Q_temp(1) = -Q_in(1);
Q_temp(2) = -Q_in(2);
Q_temp(3) = -Q_in(3);
Q_temp(4) =  Q_in(4);

V_out = Math_Qinv_Vec_Q_Mult(V_in, Q_temp);
