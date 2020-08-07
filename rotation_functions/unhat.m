function [vec] = unhat(mat)

vec(3) = mat(2,1);
vec(2) = -mat(3,1);
vec(1) = mat(3,2);

vec = vec(:);
end