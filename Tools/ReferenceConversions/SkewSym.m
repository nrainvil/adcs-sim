%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: SkewSym
% date: 10/2/2013
%
% input:        x = 3x1 vector
%
% output:       Mx_tilde = 3x3 Skew Symmetric matrix
%
% description:
% This function builds the skew symmetric matrix used for computing the
% matrix representation of the cross product.  It should be used in the
% following fashion:
% 
% A x B = SkewSym(B)'*A;
% A x B = SkewSym(A)*B;

function [Mx_tilde] = SkewSym(x)

Mx_tilde = [ 0    -x(3)  x(2);
             x(3)  0    -x(1);
            -x(2)  x(1)  0];