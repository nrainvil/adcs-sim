function [q] = dcm2q(dcm)


% compute absolute values of Euler Parameters
q1_2 = .25*(1+2*dcm(1,1)-trace(dcm));
q2_2 = .25*(1+2*dcm(2,2)-trace(dcm));
q3_2 = .25*(1+2*dcm(3,3)-trace(dcm));
q4_2 = .25*(1+trace(dcm));

[val idx] = max([q1_2 q2_2 q3_2 q4_2]);

switch idx
    case 1
        q1 = q1_2^.5;
        q2 = (dcm(1,2)+dcm(2,1))/4/q1;
        q3 = (dcm(1,3)+dcm(3,1))/4/q1;
        q4 = (dcm(2,3)-dcm(3,2))/4/q1;

    case 2
        q2 = q2_2^.5';
        q1 = (dcm(1,2)+dcm(2,1))/4/q2;
        q3 = (dcm(2,3)+dcm(3,2))/4/q2;
        q4 = (dcm(3,1)-dcm(1,3))/4/q2;
        
    case 3
        q3 = q3_2^.5;
        q1 = (dcm(1,3)+dcm(3,1))/4/q3;
        q2 = (dcm(2,3)+dcm(3,2))/4/q3;
        q4 = (dcm(1,2)-dcm(2,1))/4/q3;
        
    case 4
        q4 = q4_2^.5;
        q1 = (dcm(2,3)-dcm(3,2))/4/q4;
        q2 = (dcm(3,1)-dcm(1,3))/4/q4;
        q3 = (dcm(1,2)-dcm(2,1))/4/q4;
end

q = [q1 q2 q3 q4]';