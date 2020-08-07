% set me up cuh 

clear

% inertia matrix
J = diag([1;2;3]);
R = expm(hat(deg2rad([5;4;-6])));
J = R*J*R';

% reaction wheel jacobian
B_sc = eye(3);


% linear momentum coupling matrix
phi = [0 1 0;...
       1 0 0;
       0 .2 -.8];
 
% angular momentum coupling matrix
delta = [0 0 1;...
         0 1 0
        -.7 .1 .1];
    
% store this matrix for faster computations
T = inv(J-delta'*delta);

j = 3; % 3 modes

% damping and stiffness
zeta = [.001;.001;.001];
Delta = [.05; .2; .125] * (2*pi);

% damping and stiffness matrices 
C = zeros(j,j);
K = zeros(j,j);
for i =1:j
    C(i,i) = 2*zeta(i)*Delta(i);
    K(i,i) = Delta(i)^2;
end


           %   mrp        w                  n                       ndot 
pdot_row = [zeros(3,3) .25*eye(3)       zeros(3,j)                 zeros(3,j)];
wdot_row = [zeros(3,3) zeros(3,3)     T*delta'*K                  T*delta'*C];
ndot_row = [zeros(j,3) zeros(j,3)     zeros(j,j)                  eye(j)];
nddot_row = [zeros(j,3) zeros(j,3) (-K - delta*T*delta'*K)    (-C - delta*T*delta'*C)];

% analytical A
A_analytical = [pdot_row;wdot_row;ndot_row;nddot_row];

% analytical B
B_analytical = [zeros(3,3);
          -T*B_sc;
          zeros(j,3);
          delta*T*B_sc];

% sample time     
dt = .5;
[Ad, Bd] = c2d(A_analytical,B_analytical,dt);
[~, B_affine] = c2d(A_analytical,eye(size(Bd,1)),dt);
      

% LQR Controller 
p = 100;
d = 30;
n = 1;
ndot = 1;
r = 5;

% cost functions
Q = diag([p;p;p;d;d;d;n;n;n;ndot;ndot;ndot]);
R = diag([r;r;r]);

% LQR
K_lqr = dlqr(Ad,Bd,Q,R);

% store everything in sc struct
sc.J = J;
sc.invJ = inv(J);
sc.B = B_sc;
sc.phi = phi;
sc.delta = delta;
sc.C = C;
sc.K = K;
sc.T = T;
sc.n_modes = j;

% save
save A_and_B_3modes.mat Ad Bd K_lqr B_affine dt sc





function [A_d, B_d] = c2d(A,B,dt)
    n = size(A,1);
    p = size(B,2);

    expAB = expm([A*dt B*dt; zeros(p,n+p)]);

    A_d = expAB(1:n,1:n);
    B_d = expAB(1:n,(n+1):end);


end

function y = clamp(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
end

function [y_np1] = rk4(ODE,tn,xn,u,sc,h)
% rk4 for a single step 

xn = xn(:);

k1 = h*ODE(tn,xn,u,sc);
k2 = h*ODE(tn + h/2,xn + k1/2,u,sc);
k3 = h*ODE(tn + h/2,xn + k2/2,u,sc);
k4 = h*ODE(tn + h,xn + k3,u,sc);

y_np1 = xn + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

% MRP handling here 
p = y_np1(1:3);

% norm squared
dp2 = dot(p,p);

% if exceeds unity, swap to shadow MRP
if dp2 > 1.0
    ps = -p/dp2;
    y_np1(1:3) = ps;
end



end


function K = LQR_spacecraft(J,B_sc,p,d,r,dt)

A = [zeros(3,3) .25*eye(3);
     zeros(3,6)];
 
B = [zeros(3,3);
     -inv(J)*B_sc];
 
[~,nu] = size(B);
% dt = .5;

[Ad, Bd] = c2d(A,B,dt);

% p = 100;
% d = 30;

Q =  diag([p;p;p;d;d;d]);
R = r*eye(nu);


K = dlqr(Ad,Bd,Q,R);

end

