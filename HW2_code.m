%Connor Hughes
%4/11/21
%CH E 152B HW 1

%initialize parameters:
delta = 1;   %sample time = 1 minute
tspan = [0:delta:600];   %time span covers 600 min

%set initial condition & set point:
y_0 = 0;  
y_t = 0;

%initialize input vector:
d_1 = zeros(1, 300) + 0.5;
d_2 = zeros(1, 301) + 1;
d = cat(2, d_1, d_2);
%add measurement error:
dm = d + 0.05*randn(1, length(tspan));

%% Problem 4b: Feedforward Control
%generate feedforward control input using control law  u(t) = dm(t):
u = dm;

%simulate system:
y = zeros(1, length(tspan));
y(1) = y_0;
for i = 1:(length(tspan) - 1)
    y(i + 1) = fdfwd_tank_lvl(y(i), d(i), u(i), delta);
end

%plot tank level deviation vs time:
plot(tspan, y)
xlabel("Time (min)")
ylabel("Y (tank ht deviation)")
title("Feedforward Control with Measurement Error")

%% Problem 4c,d, and e: Feedback Control (proportional only)
%initialize proportional gain:
k_c = -0.5;

%simulate system:
y = zeros(1, length(tspan));
y(1) = y_0;
for i = 1:(length(tspan) - 1)
    y(i + 1) = fdbck_P_tank_lvl(y(i), y_t, d(i), k_c, delta);
end

%plot tank level deviation vs time:
plot(tspan, y)
xlabel("Time (min)")
ylabel("Y (tank ht deviation)")
title("Proportional Feedback Control with Measurement Error")
legend("Kc = -0.5")
%legend("Kc = -1")
%legend("Kc = -2")

%% Problem 4f: Feedback Control (PI)
k_c = -0.4;
tau_I = 1.2*delta;
int_term = zeros(1, length(tspan));

%simulate system:
y = zeros(1, length(tspan));
y(1) = y_0;
int_term(1) = y_t - y_0;
for i = 1:(length(tspan) - 1)
    y(i + 1) = fdbck_PI_tank_lvl(y(i), y_t, d(i), k_c, tau_I, int_term(i), delta);
    int_term(i + 1) = int_term(i) + delta*(y_t - y(i + 1));
end

%plot tank level deviation vs time:
plot(tspan, y)
xlabel("Time (min)")
ylabel("Y (tank ht deviation)")
title("PI Feedback Control with Measurement Error")
legend("Kc = " + k_c + ", Tau_I = " + tau_I/delta + " min")
axis([0 600 -1 1])

%% Difference Equations:
function y_kp1 = fdfwd_tank_lvl(y_k, d_k, u_k, delta)
    %compute tank level deviation at next time step:
    y_kp1 = y_k + delta*(d_k - u_k);
end

function y_kp1 = fdbck_P_tank_lvl(y_k, y_t, d_k, k_c, delta)
    %generate measured tank level deviation (with meas. error):
    y_m = y_k + 0.05*randn(1);
    %generate feedback control based on measured tank level deviation:
    u_k = k_c*(y_t - y_m);
    %compute tank level deviation at next time step:
    y_kp1 = y_k + delta*(d_k - u_k);
end

function y_kp1 = fdbck_PI_tank_lvl(y_k, y_t, d_k, k_c, tau_I, I, delta)
    %generate measured tank level deviation (with meas. error):
    y_m = y_k + 0.05*randn(1);
    %generate feedback control based on measured tank level deviation:
    u_k = k_c*(y_t - y_m) + (k_c/tau_I)*I;
    %compute tank level deviation at next time step:
    y_kp1 = y_k + delta*(d_k - u_k);
end