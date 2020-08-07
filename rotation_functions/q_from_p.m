function [q] = q_from_p(p)
% convert modified rodrigues parameter to quaternion scalar last

p = p(:);

q = (1/(1+norm(p)^2))*[2*p ;(1-norm(p)^2)];

q = q(:)/norm(q);

end
