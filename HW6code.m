%Connor Hughes
%CH E 152B HW 6

%% Exercise 14: Optimal Drone Path via MPC
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
Nx = 3;
Nu = 3;

% Constraints
% ulb = 0.3*pars.Fmax*ones(Nu,1);
% uub = pars.Fmax*ones(Nu,1);
ulb = zeros(Nu, 1);
uub = ones(Nu, 1);

xlb = [-1.75;-1.25;0];
xub = [1.75;1.25;3.3];

lb = struct('u',ulb, 'x',xlb);
ub = struct('u',uub, 'x',xub);


%%%% START: SECTION FOR STUDENTS TO MODIFY %%%%
Nt = 100; %Horizon length
Delta = 0.1; %Time step (seconds)

% Setpoint 
r_start = [-0.5;-0.5;0.5]; %x,y,z start location (meters)
r_end = [1; 0.8; 2]; %x,y,z target location (meters)

% Steady states
usp = zeros(3, 1);
%usp = pars.F_g/4*ones(Nu,1); %rotor force setpoint
xsp = zeros(Nx,1); 
xsp(1:3,:) = r_end; %target location with zero veloity


%% Insert linear discrete time model here
% Choose A,B,C,D in continuous time for velocity model
% Then ss to get a state space model;
% Then c2d on the ss ct model with Delta to obtain a discrete time model;
% then set function f = @(x,u) sys.A*(x-xsp) + sys.B*(u-usp) + xsp;

A = zeros(3);
B = eye(3);
C = eye(3);
D = zeros(3);

% P = [A, B; zeros(3), zeros(3)];
% DTsys = expm(del.*P)
% A_d = DTsys(1:3, 1:3)
% B_d = DTsys(1:3, 4:6)
% C_d = C
% D_d = D

%f = @(x,u) A_d*(x-xsp) + B_d*(u-usp) + xsp;

CTsys = ss(A, B, C, D);
del = 0.1;
DTsys = c2d(CTsys, del)

f = @(x,u) DTsys.A*(x-xsp) + DTsys.B*(u-usp) + xsp;


%% Tuning parameters
Q = eye(3);
R = eye(3);
%R = diag(1./(uub-ulb).^2);
P = Q;

% Cost functions
ell = @(x,u) (x-xsp)'*Q*(x-xsp) + (u-usp)'*R*(u-usp);
Vf = @(x) (x-xsp)'*P*(x-xsp);

% Obstacle constraint: obs(x) <= 0
obs = @(x) 0;

% Casadi functions
% replace the following ode model with
fcasadi = mpc.getCasadiFunc(f, [Nx,Nu], {'x','u'}, {'f'});
%ode_casadi = mpc.getCasadiFunc(@(x,u) crazy_ode(x,u,pars), [Nx,Nu], ...
			%{'x','u'}, {'ode'}, 'rk4', true(), 'Delta',Delta);

%%%% END: SECTION FOR STUDENTS TO MODIFY %%%%
lcasadi = mpc.getCasadiFunc(ell, {Nx, Nu}, {'x','u'}, {'l'});
Vfcasadi = mpc.getCasadiFunc(Vf, {Nx}, {'x'}, {'Vf'});
ecasadi = mpc.getCasadiFunc(obs, {Nx}, {'x'}, {'e'});

% Dimension structure and initial condition
N=struct('x',Nx,'u',Nu,'t',Nt);
x0 = zeros(Nx,1);
x0(1:3,1) = r_start;

% Create and solve nmpc problem
solver = mpc.nmpc('f',fcasadi,'l',lcasadi,'Vf',Vfcasadi,'e',ecasadi,'N',N, ...
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
plot(t,x(1,:),t,x(2,:),t,x(3,:), 'LineWidth', 1.75)
set(gca, 'FontSize', 20)
axis('tight')
legend('x','y','z')
xlabel("t (sec)")
ylabel("position (m)")

subplot(2,1,2)
stairs(t,u(1,:), 'LineWidth', 2.4)
set(gca, 'FontSize', 20)
hold on
stairs(t,u(2,:), 'LineWidth', 1.75)
stairs(t,u(3,:), 'LineWidth', 1.4)
% stairs(t,u(4,:))
axis('tight')
hold off
legend('V_x','V_y','V_z')
xlabel("t (sec)")
ylabel("velocity (m/s)")


% Save solution
%seq = [x(1:3,1:end); x(9,1:end)]';
%save -v7 "nonlineardroneMPC.mat" seq
%% Exercise 15: Identification of TCLab Hardware Model (Pulse Test)
close all; clear all; clc

% include tclab.m for initialization
tclab;

disp('Test Heater 1')
disp('LED Indicates Temperature')

figure(1)
t1s = [];
t2s = [];
h1s = [];
h2s = [];
% initial heater values
ht1 = 0;
ht2 = 0;
h1(ht1);
h2(ht2);

disp('Turn on heaters 1 and 2 to 30%')
ht1 = 30;
ht2 = ht1;
h1(ht1);
h2(ht2);
for i = 1:1000
    tic;
    t = toc;
    pause(max(0.01,1.0-t))
end
for i = 1:5000
    tic;
    if i==605
        disp('Turn up heater 1 to 40%')
        ht1 = 40;
        h1(ht1);
    end
    if i==1005
        disp('Turn both back to 30%')
        ht1 = 30;
        ht2 = 30;
        h1(ht1);
        h2(ht2);
    end
    if i==1405
        disp('Turn heater 2 to 40%')
        ht2 = 40;
        h2(ht2);
    end
    if i==1805
        disp('Turn both back to 30%')
        ht1 = 30;
        ht2 = 30;
        h1(ht1);
        h2(ht2);
    end
    if i==2205
        disp('Turn heater 2 to 40% and heater 1 to 20%')
        ht1 = 20;
        ht2 = 40;
        h1(ht1);
        h2(ht2);
    end
    if i==2605
       disp('Turn heaters 1 and 2 back to 30%')
       ht1 = 30;
       ht2 = 30;
       h1(ht1);
       h2(ht2);
    end
    if i==3005
        disp('Turn heaters 1 and 2 to 40%')
        ht1 = 40
        ht2 = 40
        h1(ht1)
        h2(ht2)
    end
    if i==3405
       disp('Turn heaters 1 and 2 back to 30%')
       ht1 = 30;
       ht2 = 30;
       h1(ht1);
       h2(ht2);
    end
    % read temperatures
    t1 = T1C();
    t2 = T2C();

    % LED brightness
    brightness1 = (t1 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness2 = (t2 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness = max(brightness1,brightness2);
    brightness = max(0,min(1,brightness)); % limit 0-1
    led(brightness);
    
    % plot heater and temperature data
    h1s = [h1s,ht1];
    h2s = [h2s,ht2];
    t1s = [t1s,t1];
    t2s = [t2s,t2];
    n = length(t1s);
    time = linspace(0,n+1,n);
    clf
    subplot(2,1,1)
    plot(time,t1s,'r.','MarkerSize',10);
    hold on
    plot(time,t2s,'b.','MarkerSize',10);
    ylabel('Temperature (degC)')
    legend('Temperature 1','Temperature 2','Location','NorthWest')
    subplot(2,1,2)
    plot(time,h1s,'r-','LineWidth',2);
    hold on
    plot(time,h2s,'b--','LineWidth',2);
    ylabel('Heater (0-5.5 V)')
    xlabel('Time (sec)')
    legend('Heater 1','Heater 2','Location','NorthWest')
    drawnow;
    t = toc;
    pause(max(0.01,1.0-t))
end

disp('Turn off heaters')
h1(0);
h2(0);

disp('Heater Test Complete')

save('TCLabID4')

%% Exercise 15 (cont.): TCLab System Identification Process
%examine and remove linear trend in the data caused by environmental disturbance
%(observed by comparing steady state values at the start and end
% of the test for the same input)
T = [t1s', t2s'];
U = [h1s', h2s'];
Tendavg = mean(T(4700:end, :));
Tstartavg = mean(T(200:600, :));
tr1 = linspace(Tstartavg(1), Tendavg(1), 5000);
tr2 = linspace(Tstartavg(2), Tendavg(2), 5000);
tr1 = tr1 - mean(tr1);
tr2 = tr2 - mean(tr2);

Tdt(:, 1) = T(:, 1) - tr1';
Tdt(:, 2) = T(:, 2) - tr2';

% Tendavg = mean(T(4700:end, :));
% Tstartavg = mean(T(200:600, :));

%plot data (with linear trend removed)
figure
subplot(2, 1, 1)
plot(time,Tdt(:, 1),'r.','MarkerSize',10);
hold on
plot(time,Tdt(:, 2),'b.','MarkerSize',10);
ylabel('Temperature (degC)')
subplot(2,1,2)
plot(time,h1s,'r-','LineWidth',2);
hold on
plot(time,h2s,'b--','LineWidth',2);
ylabel('Heater (0-5.5 V)')
xlabel('Time (sec)')
legend('Heater 1','Heater 2','Location','NorthWest')
drawnow;

%use system ID toolbox functions to estimate system from data
dat = iddata(Tdt(600:end, :), U(600:end, :), 1.0004);
tsim = time(600:end);
G = tfest(dat, 1);
[pv, pvsd] = getpvec(G)
G2 = tfest(dat, 2);
[pv2, pvsd2] = getpvec(G2);
Gss = ssest(dat, 5);
Gssd = c2d(Gss, 1.0004);

%simulate solution of estimated system for the same input used in the pulse
%test and plot on the same axes as the test data, for comparison
subplot(2, 1, 1)
x0 = Gss.C\(Tdt(600, :)');
y = lsim(Gssd, U(600:end, :), tsim, x0, 'zoh');
plot(tsim, y, 'g', 'LineWidth', 1)
legend('Temp. 1 Detrended','Temp. 2 Detrended', 'Temp. 1 Simulated', 'Temp. 2 Simulated', 'Location','NorthWest')
save('TCLabID4comp')

%% Exercise 16: PID Control with Estimated TCLab Model
%close all;
load('TCLabID4Comp')

%Find steady-state gains for estimated 2x2 TF matrix
K_11 = G.Numerator{1, 1}/G.Denominator{1, 1}(2);
K_12 = G.Numerator{1, 2}/G.Denominator{1, 2}(2);
K_21 = G.Numerator{2, 1}/G.Denominator{2, 1}(2);
K_22 = G.Numerator{2, 2}/G.Denominator{2, 2}(2);

K = [K_11, K_12; K_21, K_22];
RGA = K.*(inv(K)');

%Compute time constants for each TF:
tau_11 = 1/G.Denominator{1, 1}(2);
tau_12 = 1/G.Denominator{1, 2}(2);
tau_21 = 1/G.Denominator{2, 1}(2);
tau_22 = 1/G.Denominator{2, 2}(2);

%PID Tuning:
K_c1 = 1/K_11;
K_c2 = 1/K_22;
Tau_I1 = tau_11;
Tau_I2 = tau_22;

K_c1 = 1.4*K_c1;
K_c2 = 1.3*K_c2;
Tau_I1 = 0.65*Tau_I1;
Tau_I2 = 1*Tau_I2;

tsam = 1;
dsys = c2d(ss(G), tsam);
ad = dsys.a;
bd = dsys.b;
cd = dsys.c;

%generate set points for closed-loop system:
nsim = 18000;
T1_sp = zeros(1, 18000) + 40;
T1_sp(14000:15000) = 45;
T1_sp(15000:16000) = 40;
T1_sp(16000:17000) = 45;
T1_sp(17000:end) = 40;
T2_sp = zeros(1, 18000) + 40;
T2_sp(14000:15000) = 50;
T2_sp(15000:16000) = 45;
T2_sp(16000:17000) = 40;
T2_sp(17000:end) = 40;

%compute system response to above set points:
T_sp = [T1_sp; T2_sp]';
y0 = [40; 40];
x0 = dsys.C\y0;
time=(0:nsim-1)*tsam;
order = length(ad(:,1));
x = zeros(order,nsim);
x(1, :) = x(1, :) + x0(1);
x(2, :) = x(2, :) + x0(2);
x(3, :) = x(3, :) + x0(3);
x(4, :) = x(4,  :) + x0(4);
interr = zeros(2,nsim);
trackerr = zeros(2,nsim);
y = zeros(2, nsim);
y(1, :) = y(1, :) + y0(1);
y(2, :) = y(2, :) + y0(2);
y_t = [T1_sp; T2_sp];
u(1, 1) = 36;
u(2, 1) = 24;

for k = 2:nsim
    y(:,k) = cd*x(:,k) + 0.0*randn(2,1);
    trackerr(:,k) = y_t(:,k)-y(:,k);
    u(1,k) = K_c1*(trackerr(1,k)+ 1/Tau_I1*interr(1,k));
    u(2,k) = K_c2*(trackerr(2,k)+ 1/Tau_I2*interr(2,k));
    if(k == nsim) 
        break 
    end
    interr(:,k+1) = interr(:,k)+trackerr(:,k)*tsam;
    x(:,k+1) = ad*x(:,k)+bd*u(:,k);
end

figure()
subplot(2,1,1)
ylabel("y", "rotation", 0)
plot(time, y, time, y_t)
%plot(time(14000:end), y(:, 1400:end), time(14000:end), y_t(:, 14000:end))
ylim([0 60])
xlabel('time (s)')
ylabel('temperature (deg C)')
legend('Temperature 1', 'Temperature 2', 'Temp. 1 Set Point', 'Temp. 2 Set Point')
title("P Gains: " + PI_1.Kp + " " + PI_2.Kp + ",  I Gains: " + PI_1.Ki + " " + PI_2.Ki)
%title("Temperatures and Set Points");

subplot(2,1,2)
ylabel ("u", "rotation", 0)
stairs(time, u')
xlabel('time (s)')
ylabel('Heater Setting (% Full Voltage)')
legend('Heater 1', 'Heater 2')
title("P Gains: " + K_c1 + " " + K_c2 + ",  Taus: " + Tau_I1 + " " + Tau_I2)
%title("Heater Settings")

% ax = gca
% ax.FontSize = 20;

%% Step tests 
step(Gssd)
step(G)

%% Test

%PID Tuning:
K_c1 = 1/K_11;
K_c2 = 1/K_22;
Tau_I1 = tau_11;
Tau_I2 = tau_22;

PI_1 = pid(K_c1, 1*K_c1/Tau_I1);
PI_2 = pid(K_c2, 1*K_c2/Tau_I2);

C = [PI_1, 0; 0, PI_2];

GC = G*C;

tsam = 1;

sys = feedback(GC, eye(2));
dsys = c2d(ss(sys), tsam);

%generate set points for closed-loop system:
time = linspace(0, 7999, 8000);

T1_sp = zeros(1, 8000) + 40;

T2_sp = zeros(1, 8000) + 40;

%plot set points and system response
figure
y_sp = [T1_sp; T2_sp]';
y = lsim(dsys, y_sp, time);
plot(time, y)
hold on
plot(time, y_sp)
ylim([0 60])
xlabel('time (s)')
ylabel('temperature (deg C)')
legend('Temperature 1', 'Temperature 2', 'Temp. 1 Set Point', 'Temp. 2 Set Point')
title("P Gains: " + PI_1.Kp + " " + PI_2.Kp + ",  I Gains: " + PI_1.Ki + " " + PI_2.Ki)

e = y_sp - y;
u = lsim(tst, e, time);


figure
plot(time, u)
% hold on
% plot(time, y_sp)
ylim([0 60])
xlabel('time (s)')
ylabel('Heater Setting (% Full Voltage)')
legend('Heater 1', 'Heater 2')
title("P Gains: " + PI_1.Kp + " " + PI_2.Kp + ",  I Gains: " + PI_1.Ki + " " + PI_2.Ki)

%% JBR-code
    clear u
    tsam = 1;  
    ad = dsys.a;
    bd = dsys.b;
    cd = dsys.c;

    nsim = 8000;
    T1_sp = zeros(1, 8000) + 40;
    T1_sp(4000:5000) = 45;
    T1_sp(5000:6000) = 40;
    T1_sp(6000:7000) = 45;
    T1_sp(7000:end) = 40;

    T2_sp = zeros(1, 8000) + 40;
    T2_sp(4000:5000) = 50;
    T2_sp(5000:6000) = 45;
    T2_sp(6000:7000) = 40;
    T2_sp(7000:end) = 40;

    time=(0:nsim-1)*tsam;
    order = length(ad(:,1));
    x = zeros(order,nsim);
    interr = zeros(2,nsim);
    y = zeros(2, nsim);
    y_t = [T1_sp; T2_sp];
    y_t(:,1) = 0;

    for k = 2:nsim
        y(:,k) = cd*x(:,k) + 0.0*randn(2,1);
        trackerr(:,k) = y_t(:,k)-y(:,k);
        u(1,k) = K_c1*( trackerr(1,k)+ 1/Tau_I1*interr(1,k) );
        u(2,k) = K_c2*( trackerr(2,k)+ 1/Tau_I2*interr(2,k) );
        if(k == nsim) 
            break 
        end
        interr(:,k+1) = interr(:,k)+trackerr(:,k)*tsam;
        x(:,k+1) = ad*x(:,k)+bd*u(:,k);
    end
    figure()
    subplot(2,1,1)
    ylabel("y", "rotation", 0)
    plot(time, y, time, y_t)
    xlabel('time (s)')
    ylabel('temperature (deg C)')
    legend('Temperature 1', 'Temperature 2', 'Temp. 1 Set Point', 'Temp. 2 Set Point')
    title("P Gains: " + PI_1.Kp + " " + PI_2.Kp + ",  I Gains: " + PI_1.Ki + " " + PI_2.Ki)
    
    subplot(2,1,2)
    ylabel ("u", "rotation", 0)
    stairs(time, u')
    xlabel('time (s)')
    ylabel('Heater Setting (% Full Voltage)')
    legend('Heater 1', 'Heater 2')
    title("P Gains: " + PI_1.Kp + " " + PI_2.Kp + ",  I Gains: " + PI_1.Ki + " " + PI_2.Ki)

