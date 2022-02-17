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
Nt = 100; %Horizon length
Delta = 0.1; %Time step (seconds)

% Setpoint 
r_start = [-0.5;-0.5;0.5]; %x,y,z start location (meters)
r_end = [1; 0.8; 2]; %x,y,z target location (meters)

% Steady states
usp = pars.F_g/4*ones(Nu,1); %rotor force setpoint
xsp = zeros(Nx,1); 
xsp(1:3,:) = r_end; %target location with zero veloity


%% Insert linear discrete time model here
%% Choose A,B,C,D in continuous time for velocity model
%% Then ss to get a state space model;
%% Then c2d on the ss ct model with Delta to obtain a discrete time model;
%% then set function f = @(x,u) sys.A*(x-xsp) + sys.B*(u-usp) + xsp;
%%

%% Tuning parameters
Q = diag([ones(6,1);zeros(6,1)]);
R = diag(1./(uub-ulb).^2);
P = Q;

% Cost functions
ell = @(x,u) (x-xsp)'*Q*(x-xsp) + (u-usp)'*R*(u-usp);
Vf = @(x) (x-xsp)'*P*(x-xsp);

% Obstacle constraint: obs(x) <= 0
obs = @(x) 0;

% Casadi functions
%% replace the following ode model with
%% fcasadi = mpc.getCasadiFunc(f, [Nx,Nu], {'x','u'}, {'f'});
ode_casadi = mpc.getCasadiFunc(@(x,u) crazy_ode(x,u,pars), [Nx,Nu], ...
			{'x','u'}, {'ode'}, 'rk4', true(), 'Delta',Delta);

%%%% END: SECTION FOR STUDENTS TO MODIFY %%%%
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
figure(1)
subplot(2,1,1)
plot(t,x(1,:),t,x(2,:),t,x(3,:))
axis('tight')
legend('r_1 (x)','r_2 (y)','r_3 (z)')
xlabel("t (sec)")
ylabel("r (m)")

subplot(2,1,2)
stairs(t,u(1,:))
hold on
stairs(t,u(2,:))
stairs(t,u(3,:))
stairs(t,u(4,:))
axis('tight')
hold off
legend('F_1','F_2','F_3','F_4')
xlabel("t (sec)")
ylabel("F (N)")


% Save solution
%seq = [x(1:3,1:end); x(9,1:end)]';
%save -v7 "nonlineardroneMPC.mat" seq