function [phi] = phi_from_q(quat)
quat = quat(:);

theta = 2*atan2(norm(quat(1:3)),quat(4));
if theta>pi
    theta = theta-2*pi;
end

r = (quat(1:3))/norm(quat(1:3));

phi = r*theta;
phi = phi(:);

end