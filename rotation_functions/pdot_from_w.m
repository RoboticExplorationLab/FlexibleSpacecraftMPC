function [pdot] = pdot_from_w(p,w)
% this is the kinematics of the modified rodrigues parameter assuming that
% attitude is being denoted as N_R_B using the kane/levinson convention
p = p(:);
w = w(:);

pdot = ((1+norm(p)^2)/4)*(eye(3) + 2*(hat(p)^2 + hat(p))/(1+norm(p)^2))*w;

end
