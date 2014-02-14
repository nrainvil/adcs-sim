function [ out_vec ] = ortho_norm( in_vec1, in_vec2 )
%ORTHO_NORM
    in_vec1_h = in_vec1/norm(in_vec1);
    in_vec2_h = in_vec2/norm(in_vec2);
    out_vec = in_vec1_h - dot(in_vec1_h,in_vec2_h)*in_vec2_h;
    out_vec = out_vec/norm(out_vec);
end

