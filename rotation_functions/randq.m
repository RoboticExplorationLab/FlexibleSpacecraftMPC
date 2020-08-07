function [q] = randq

q = rand(4,1);
q = q/norm(q);

end
