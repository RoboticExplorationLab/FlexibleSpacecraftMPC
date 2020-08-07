%% run monte carlo runs for uncertainty about modes 
clear


N = 10;

mc_runs_mpc = cell(N,1);
mc_runs_lqr = cell(N,1);
MSE_mpc = zeros(N,1);
MSE_lqr = zeros(N,1);


for i = 1:N

[t_vec,angle_error_mpc_sim,angle_error_lqr] = run_mc();
mc_runs_mpc{i} = angle_error_mpc_sim;
mc_runs_lqr{i} = angle_error_lqr;

MSE_mpc(i) = mean(angle_error_mpc_sim.^2);
MSE_lqr(i) = mean(angle_error_lqr.^2);
% plot(t_vec(1:end-1),rad2deg(angle_error_mpc_sim))
i
end

%% plotting



mpc_mat = cell2mat(mc_runs_mpc);
lqr_mat = cell2mat(mc_runs_lqr);


mpc_lower = zeros(200,1);
mpc_upper = zeros(200,1);
lqr_upper = zeros(200,1);
lqr_lower = zeros(200,1);
mpc_avg = mpc_lower;
mpc_std = mpc_avg;
lqr_std = mpc_std;
lqr_avg = mpc_avg;

for i = 1:200
    mpc_lower(i) = min(mpc_mat(:,i));
    mpc_upper(i) = max(mpc_mat(:,i));
    lqr_lower(i) = min(lqr_mat(:,i));
    lqr_upper(i) = max(lqr_mat(:,i));
end

for i = 1:200
    mpc_avg(i) = min(mpc_mat(:,i));
    mpc_std(i) = std(mpc_mat(:,i));
    lqr_avg(i) = min(lqr_mat(:,i));
    lqr_std(i) = std(lqr_mat(:,i));
end


col1 = [.25 .25 .25];
col2 = [.8 .8 .8];

x = t_vec(1:end-1);


figure
hold on 
x2 = [x,fliplr(x)];
inBetween_mpc = [(mpc_avg - 3*mpc_std)', fliplr((mpc_avg + 3*mpc_std)')];
inBetween_lqr = [(lqr_avg - 3*lqr_std)', fliplr((lqr_avg + 3*lqr_std)')];

h(2) = fill(x2,rad2deg(inBetween_lqr),col2)
plot(x,rad2deg((lqr_avg - 3*lqr_std)),'Color',col2)% comment this out to switch
plot(x,rad2deg((lqr_avg + 3*lqr_std)),'Color',col2)% comment this out to switch
h(1) = fill(x2,rad2deg(inBetween_mpc),col1)
plot(x,rad2deg((mpc_avg - 3*mpc_std)),'Color',col1)
plot(x,rad2deg((mpc_avg + 3*mpc_std)),'Color',col1)
% plot(x,rad2deg((lqr_avg - 3*lqr_std)),'Color',col2) % comment htis out to
% plot(x,rad2deg((lqr_avg + 3*lqr_std)),'Color',col2) % comment this out
ylim([0 10])
title('Monte-Carlo 3-\sigma Pointing Error Performance (1000 Trials)')
% legend('MPC','LQR')
[~,hObj]=legend(h, 'MPC','LQR');  
ylabel('Pointing Error (deg)')
xlabel('Time (s)')
hold off 

%%

figure
hold on 

col1 = [.25 .25 .25];
col2 = [.7 .7 .7];

for i = 1:N
    plot(t_vec(1:end-1),round(rad2deg(mc_runs_lqr{i}),2),'Color',col2)
end

for i = 1:N
    plot(t_vec(1:end-1),round(rad2deg(mc_runs_mpc{i}),2),'Color',col1)
end

h(1) = plot(NaN,NaN,'Color',col1);
h(2) = plot(NaN,NaN,'Color',col2);
%legend(h, 'Run 1','Run 100');
[~,hObj]=legend(h, 'MPC','LQR');           % return the handles array
hL=findobj(hObj,'type','line');  % get the lines, not text
set(hL,'linewidth',6) 

title('Model Uncertainty Monte-Carlo Pointing Error (1000 Trials)')
ylabel('Pointing Error (degrees)')
xlabel('Time (seconds)')
hold off

%%

hist_mpc = rad2deg(sqrt(MSE_mpc));
hist_lqr = rad2deg(sqrt(MSE_lqr));

figure

hold on 

title('Root-Mean-Square Pointing Error Monte-Carlo (1000 Trials)')
histogram(hist_mpc,20,'FaceColor',[.1 .1 .1])
histogram(hist_lqr,20,'FaceColor','w')
xlim([0 3])
% ylim([0 220])
ylabel('Frequency')
xlabel('Root-Mean-Square Pointing Error (degrees)')

legend('MPC','LQR')

hold off 

%% max stuff 

max_mpc = zeros(N,1);
max_lqr = zeros(N,1);

for i = 1:N

max_mpc(i) = max(mc_runs_mpc{i});
max_lqr(i) = max(mc_runs_lqr{i});

end


figure
hold on 
title('Maximum Pointing Error Monte-Carlo (1000 Trials)')
histogram(rad2deg(max_mpc),15,'FaceColor',[.1 .1 .1])
histogram(rad2deg(max_lqr),15,'FaceColor','w')
xlim([0 12])
% ylim([0 210])
ylabel('Frequency')
xlabel('Maximum Pointing Error (Degrees)')
legend('MPC','LQR')

%% supporting functions 
function [t_vec,angle_error_mpc_sim,angle_error_lqr] = run_mc()

[Ad, Bd, K_lqr, B_affine, dt, sc] = setup_linearized_model();

[t_vec,angle_error_mpc_sim,angle_error_lqr] = run_sim(Ad, Bd, K_lqr, B_affine, dt, sc);

end

function [Ad, Bd, K_lqr, B_affine, dt, sc] = setup_linearized_model()

J = diag([1;2;3]);
R = expm(hat(deg2rad([5;4;-6])));
J = R*J*R';

B_sc = eye(3);



phi = normalize([0;1;0],'norm')';
delta = normalize([0;0;1],'norm')';
T = inv(J-delta'*delta);

zeta = .001;
Delta = .05 * (2*pi);
C = 2*zeta*Delta;
K = Delta^2;

j = size(C,1); % modal coordiante number 
           %   mrp        w          n         ndot 
pdot_row = [zeros(3,3) .25*eye(3) zeros(3,j) zeros(3,j)];
wdot_row = [zeros(3,3) zeros(3,3)   T*delta'*K    T*delta'*C];
ndot_row = [zeros(j,3) zeros(j,3) zeros(j,j) eye(j)];
nddot_row = [zeros(j,3) zeros(j,3) (-K - delta*T*delta'*K)    (-C - delta*T*delta'*C)];

A_anal = [pdot_row;wdot_row;ndot_row;nddot_row];

B_anal = [zeros(3,3);
          -T*B_sc;
          zeros(j,3);
          delta*T*B_sc];

      
dt = .5;
[Ad, Bd] = c2d(A_anal,B_anal,dt);
[~, B_affine] = c2d(A_anal,eye(8),dt);
      

% LQR Controller 
p = 100;
d = 30;
r = 5;

K_lqr = LQR_spacecraft(J,B_sc,p,d,r,dt);

sc.J = J;
sc.invJ = inv(J);
sc.B = B_sc;
sc.phi = phi;
sc.delta = delta;
sc.C = C;
sc.K = K;
sc.T = T;


% fake stuff for the sim s
zeta = .001*(1+.1*randn);
Delta = .05 * (2*pi)*(1+.1*randn);
C_fake = 2*zeta*Delta;
K_fake = Delta^2;

sc.C_fake = C_fake;
sc.K_fake = K_fake;

err1 = deg2rad(5)*randn*normalize(randn(3,1),'norm');
err2 = deg2rad(5)*randn*normalize(randn(3,1),'norm');
sc.phi_fake = (1+.05*randn)*phi*expm(hat(err1));
sc.delta_fake = (1+.05*randn)*delta*expm(hat(err2));

% save A_and_B_shifted_up.mat Ad Bd K_lqr B_affine dt sc

end


function y = clamp(x,bl,bu)
  % return bounded value clipped between bl and bu
  y=min(max(x,bl),bu);
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


function [t_vec,angle_error_mpc_sim,angle_error_lqr] = run_sim(Ad, Bd, K_lqr, B_affine, dt, sc)
[nx, nu] = size(Bd);

% Constraints
umin = -.01*ones(3,1);
umax = -umin;
xmin = -Inf(8,1);
xmax = -xmin;

w.mrp = 20;
w.w = 1;
w.n = 1;
w.nd = 1;

% Objective function
Q = diag([w.mrp*ones(3,1);w.w*ones(3,1);w.n;w.nd]);
QN = Q;
R = 0.0*eye(nu);

% Initial and reference states
x0 = [0;0;.0;zeros(5,1)];
xr = zeros(8,1);

% Prediction horizon
N = 50;

% Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
% - quadratic objective
P = blkdiag( kron(speye(N), Q), QN, kron(speye(N), R) );
% - linear objective
q = [repmat(-Q*xr, N, 1); -QN*xr; zeros(N*nu, 1)];
% - linear dynamics
Ax = kron(speye(N+1), -speye(nx)) + kron(sparse(diag(ones(N, 1), -1)), Ad);
Bu = kron([sparse(1, N); speye(N)], Bd);
Aeq = [Ax, Bu];

% add the affine stuff here
sc.tau= 10*[.01;.01;.01];
sc.a = [.1;.1;.1];

% throw a disturbance in there at 7.5 seconds 
affine_matrix = zeros(nx,N);
affine_instance = affine_fx(sc.tau,sc.a,sc);
affine_matrix(:,15) = B_affine*affine_instance;


leq = [-x0; -vec(affine_matrix)];
ueq = leq;
% - input and state constraints
Aineq = speye((N+1)*nx + N*nu);
lineq = [repmat(xmin, N+1, 1); repmat(umin, N, 1)];
uineq = [repmat(xmax, N+1, 1); repmat(umax, N, 1)];
% - OSQP constraints
A = [Aeq; Aineq];
l = [leq; lineq];
u = [ueq; uineq];

% Create an OSQP object
prob = osqp;

% Setup workspace
prob.setup(P, q, A, l, u,...
           'eps_abs',1e-6,...
           'eps_rel',1e-6,...
           'eps_prim_inf',1e-6,...
           'eps_dual_inf', 1e-6,...
           'verbose',false);

% solve 
res = prob.solve();

%% plotting 

% Simulate in closed loop with true nonlinear dynamics
N_sim = 200;
% dt = .5;
t_vec = 0:dt:N_sim*dt;
ut_vec = t_vec(1:end-1);
X_mpc_sim = zeros(11,N_sim);
X_mpc_sim(:,1) = [x0;0;0;0];
U_mpc_sim = zeros(nu,N_sim-1);
X_lqr_sim = zeros(11,N_sim);
U_lqr_sim = zeros(nu,N_sim-1);


for i = 1 : N_sim
    %% MPC simulation
    % Solve the QP
    res = prob.solve();

    % Check solver status
%     if ~strcmp(res.info.status, 'solved')
%         error('OSQP did not solve the problem!')
%     end

    % Apply first control input to the plant
    U_mpc_sim(:,i) = clamp(res.x((N+1)*nx+1:(N+1)*nx+nu),umin(1),umax(1));

    % step forward dt
    X_mpc_sim(:,i+1) = rk4(@ODE,t_vec(i),X_mpc_sim(:,i),U_mpc_sim(:,i),sc,dt);

    % Update initial state
    [l,u] = shift_affine(X_mpc_sim(1:8,i+1),l,u,N, nx);
    prob.update('l', l, 'u', u);
    
    %% LQR simulation
    
    % control input 
    U_lqr_sim(:,i) = clamp(-K_lqr*X_lqr_sim(1:6,i),umin(1),umax(1));
    
    % step forward dt
    X_lqr_sim(:,i+1) = rk4(@ODE,t_vec(i),X_lqr_sim(:,i),U_lqr_sim(:,i),sc,dt);
    
    
end

%%

angle_error_mpc_sim = zeros(1,N_sim);
angle_error_lqr = zeros(1,N_sim);
for i = 1:N_sim
     angle_error_lqr(i) = norm(phi_from_p(X_lqr_sim(1:3,i)));
    angle_error_mpc_sim(i) = norm(phi_from_p(X_mpc_sim(1:3,i)));
end


%% Plot both approaches 
% 
% angle_error_mpc = angle_error_mpc_sim;
% X = X_mpc_sim;
% U = U_mpc_sim;
% figure
% hold on 
% sgtitle('MPC vs LQR Station-Keeping')
% % mrp 
% % subplot(3,2,1)
% % plot(t_vec,X(1:3,:)')
% % title('MPC MRP')
% % legend('p_1','p_2','p_3')
% % ylabel('MRP')
% % xlabel('Time (s)')
% % ylim([-.02 .08])
% % hold off 
% % 
% % subplot(3,2,2)
% % plot(t_vec,X_sim(1:3,:)')
% % title('LQR MRP')
% % legend('p_1','p_2','p_3')
% % ylabel('MRP')
% % xlabel('Time (s)')
% % ylim([-.02 .08])
% % hold off 
% subplot(3,2,1)
% plot(ut_vec,rad2deg(angle_error_mpc))
% title('MPC Pointing Error')
% ylabel('Pointing Error (deg)')
% xlabel('Time (s)')
% ylim([0 10])
% hold off 
% 
% subplot(3,2,2)
% plot(ut_vec,rad2deg(angle_error_lqr))
% title('LQR Pointing Error')
% ylabel('Pointing Error (deg)')
% xlabel('Time (s)')
% ylim([0 10])
% hold off 
% 
% % modal coordinate 
% subplot(3,2,3)
% plot(t_vec,X_mpc_sim(7,:)')
% title('MPC Modal Coordinate')
% ylabel('Modal Coordinate')
% xlabel('Time (s)')
% ylim([-.2 .2])
% hold off 
% 
% subplot(3,2,4)
% plot(t_vec,X_lqr_sim(7,:)')
% title('LQR Modal Coordinate')
% ylabel('Modal Coordinate')
% xlabel('Time (s)')
% ylim([-.2 .2])
% hold off 
% 
% % Control 
% subplot(3,2,5)
% stairs(ut_vec,U_mpc_sim')
% title('MPC U')
% ylabel('N*m')
% xlabel('Time (s)')
% legend('u_1','u_2','u_3')
% ylim([-.011 .011])
% hold off 
% 
% subplot(3,2,6)
% stairs(ut_vec,U_lqr_sim(1:3,:)')
% title('LQR U')
% ylabel('N*m')
% xlabel('Time (s)')
% legend('u_1','u_2','u_3')
% ylim([-.011 .011])
% hold off 
% set(gcf,'Position',[100 100 1700 900])

end


% function y = clamp(x,bl,bu)
%   % return bounded value clipped between bl and bu
%   y=min(max(x,bl),bu);
% end

function aff = affine_fx(tau,a,sc)

delta = sc.delta;
phi = sc.phi;
T = sc.T;

aff = [zeros(3,1);
                 (T*tau + T*delta'*phi*a);
                 zeros(1);
                (-delta*T*tau - delta*T*delta'*phi*a - phi*a)];
end


function xdot = ODE(t,x,u,sc)

p    = x(1:3);
w    = x(4:6);
n    = x(7);
ndot = x(8);
r    = x(9:11);

% tau = u.tau;
% a = u.a; 
% rdot = u.rdot;
if t>= 7.0 && t <7.5
        tau = sc.tau;
        a   = sc.a;
else
        tau = [0.0;0.0;0.0];
        a   = [0.0;0.0;0.0];
end

wdot = sc.T*(tau - sc.B*u  - ...
        cross(w,sc.J*w + sc.delta_fake'*ndot + sc.B*r) +...
        sc.delta_fake'*(sc.C_fake*ndot + sc.K_fake*n + sc.phi_fake*a));

pdot = pdot_from_w(p,w);

nddot = -sc.delta_fake*wdot -sc.C_fake*ndot - sc.K_fake*n - sc.phi_fake*a;

xdot = [pdot;wdot;ndot;nddot;u];


end





function [A_d, B_d] = c2d(A,B,dt)
    n = size(A,1);
    p = size(B,2);

    expAB = expm([A*dt B*dt; zeros(p,n+p)]);

    A_d = expAB(1:n,1:n);
    B_d = expAB(1:n,(n+1):end);


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

% if MRP exceeds unity, swap to shadow MRP
if dp2 > 1.0
    ps = -p/dp2;
    y_np1(1:3) = ps;
end



end

function [l,u] = shift_affine(x0,l,u,N, nx)

leq   = l(1:(N+1)*nx);
lineq = l((N+1)*nx+1:end);
% ueq   = u(1:(N+1)*nx);
uineq = u((N+1)*nx+1:end);

leq_shifted = zeros(size(leq));

for i = 1:length(leq_shifted)-nx
    leq_shifted(i) = leq(i+nx);
end

leq_shifted(1:nx) = -x0;


ueq_shifted = leq_shifted;


l = [leq_shifted; lineq];
u = [ueq_shifted; uineq];


end