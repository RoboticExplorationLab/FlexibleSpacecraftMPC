function phi = phi_from_p(p)

if norm(p) == 0 
    phi = zeros(3,1);
else
    q = q_from_p(p);
    phi = phi_from_q(q);
end


end
