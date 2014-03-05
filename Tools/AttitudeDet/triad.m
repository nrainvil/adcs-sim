function [ R ] = triad( s_b,s_i,b_b,b_i )
%TRIAD - Generate Deterministic Rotation Matrix 
% from sun and mag vectors
% input: s_b - Sun Body Vector, s_i Sun Inertial Vector
%       b_b - Mag Body Vector, b_i Mag Inertial Vector
% output: R - 3x3 Cosine matrix

    t1_b = s_b;
    t1_i = s_i;
    t2_b = cross(s_b,b_b)/norm(cross(s_b,b_b));
    t2_i = cross(s_i,b_i)/norm(cross(s_i,b_i));
    t3_b = cross(t1_b,t2_b);
    t3_i = cross(t1_i,t2_i);
    
    R = [t1_b,t2_b,t3_b]*[t1_i,t2_i,t3_i]';
end

