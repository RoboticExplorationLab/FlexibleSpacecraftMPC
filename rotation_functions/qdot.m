function [product] = qdot(q1,q2)
%scalar last quaternion multiplication, hamilton product, qdot

%make column vectors
q1=q1(:);
q2 = q2(:);

%break quats into s and v 
v1 = q1(1:3);
s1 = q1(4);
v2 = q2(1:3);
s2 = q2(4);

%multiply 
product = zeros(4,1);
product(1:3) = s1*v2 + s2*v1 + cross(v1,v2);
product(4) = s1*s2 - v1'*v2;
product = product(:);
end