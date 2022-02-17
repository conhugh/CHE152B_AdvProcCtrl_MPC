mpc = import_mpctools(); %Import mpctools

% Crazyflie Parameters.
pars = struct();
pars.m = 0.033;            % kg      mass
pars.F_g =  pars.m * 9.81; % N       weight
pars.b = 0.046;            % m       half width of drone
pars.J_xx = 1.9e-5;        % kg*m^2  moment of inertia about x-axis
pars.J_yy = 1.9e-5;        % kg*m^2  moment of inertia about y-axis
pars.J_zz = 2.6e-5;        % kg*m^2  moment of inertia about z-axis
pars.d_x =  0.005;         % N/(m/s) drag coeff. (translat. motion)
pars.d_y =  0.005;         % N/(m/s) drag coeff.
pars.d_z =  0.01;          % N/(m/s) drag coeff.
pars.kappa = 0.0037;       % m       ratio k_T/k_F (torque and force coeff.)
pars.Fmax = 0.12;	   % N       maximum force of a rotor

% Dimensions for the nonlinear system
Nx = 12;
Nu = 4;

% Constraints
ulb = 0.3*pars.Fmax*ones(Nu,1);
uub = pars.Fmax*ones(Nu,1);

xlb = [-1.75;-1.25;0;-1;-1;-1;-0.2;-0.2;-pi;-pi/2;-pi/2;-pi/2];
xub = [1.75;1.25;3.3; 1; 1; 1; 0.2; 0.2; pi; pi/2; pi/2; pi/2];

lb = struct('u',ulb, 'x',xlb);
ub = struct('u',uub, 'x',xub);


%%%% START: SECTION FOR STUDENTS TO MODIFY %%%%
Nt = 80; %Horizon length
Delta = 0.1; %Time step (seconds)

% Setpoint 
r_start = [-1;0;0.5]; %x,y,z start location (meters)
r_end = [1; 0; 0.5]; %x,y,z target location (meters)

% Steady states
usp = pars.F_g/4*ones(Nu,1); %rotor force setpoint
xsp = zeros(Nx,1); 
xsp(1:3,:) = r_end; %target location with zero velocity

% Tuning parameters
Qscl = 5;
Rscl = 1;
Q = diag([Qscl*ones(6,1);zeros(6,1)]);
R = diag(Rscl./(uub-ulb).^2);
P = Q;

% Cost functions
ell = @(x,u) (x-xsp)'*Q*(x-xsp) + (u-usp)'*R*(u-usp);
Vf = @(x) (x-xsp)'*P*(x-xsp);

% Obstacle constraint: obs(x) <= 0
%obs = @(x) 0;
obs = @(x) (-10*x(1)^2 + 1.7) - x(3);

%%%% END: SECTION FOR STUDENTS TO MODIFY %%%%

% Casadi functions
ode_casadi = mpc.getCasadiFunc(@(x,u) crazy_ode(x,u,pars), [Nx,Nu], ...
			{'x','u'}, {'ode'}, 'rk4', true(), 'Delta',Delta);
lcasadi = mpc.getCasadiFunc(ell, {Nx, Nu}, {'x','u'}, {'l'});
Vfcasadi = mpc.getCasadiFunc(Vf, {Nx}, {'x'}, {'Vf'});
ecasadi = mpc.getCasadiFunc(obs, {Nx}, {'x'}, {'e'});

% Dimension structure and initial condition
N=struct('x',Nx,'u',Nu,'t',Nt);
x0 = zeros(Nx,1);
x0(1:3,1) = r_start;

% Create and solve nmpc problem
solver = mpc.nmpc('f',ode_casadi,'l',lcasadi,'Vf',Vfcasadi,'e',ecasadi,'N',N, ...
			'lb',lb,'ub',ub,'x0',x0,'verbosity',3);
solver.solve();

% Extract optimal solution
x = solver.var.x;
u = solver.var.u;
u = [u,u(:,end)];
t = (0:Nt)*Delta;

%Plot
figure
subplot(2,2,1)
plot(t,x(1,:),t,x(2,:),t,x(3,:), 'linewidth', 1.2)
ylim([-1.2, 1.8])
legend('r_1 (x)','r_2 (y)','r_3 (z)', 'FontSize', 12)
xlabel("t (sec)")
ylabel("position (m)")
title('Position Trajectory')
ax = gca
ax.FontSize = 16;

subplot(2,2,2)
plot(t,x(4,:),t,x(5,:),t,x(6,:),'linewidth', 1.2)
ylim([-1.0, 1.1])
legend('v_x','v_y','v_z', 'FontSize', 12)
xlabel("t (sec)")
ylabel("velocity (m/s)")
title('Velocity Trajectory')
ax = gca
ax.FontSize = 16;

subplot(2,2,3)
stairs(t,u(1,:), 'linewidth', 1.2)
hold on
stairs(t,u(2,:), 'linewidth', 1.2)
stairs(t,u(3,:), 'linewidth', 1.2)
stairs(t,u(4,:), 'linewidth', 1.2)
axis('tight')
ylim([0.02, 0.13])
hold off
legend('F_1','F_2','F_3','F_4', 'FontSize', 12)
xlabel("t (sec)")
ylabel("Force (N)")
title('Rotor Thrust Trajectories')

ax = gca
ax.FontSize = 16;

ob = @(x) -10*x(1)^2 + 1.5;
obstacle = zeros(1, length(x(1, :)));
for i = 1:length(x(1, :))
    obstacle(i) = ob(x(:, i));
end

subplot(2, 2, 4)
plot(x(1, :),x(3,:), x(1, :), obstacle)
ylim([0 2]);
legend('Drone Flight Path','Obstacle', 'FontSize', 12)
xlabel("x (m)")
ylabel("z (m)")
title('Path Visualization')
ax = gca
ax.FontSize = 16;

%Save solution
seq = [x(1:3,1:end); x(9,1:end)]';
save("ConnorHughesDroneTrajec.mat", 'seq')