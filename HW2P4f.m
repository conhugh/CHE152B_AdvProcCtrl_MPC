%% Problem 4f: Feedback Control (PI)
k_c = -0.4;
tau_I = 1.2*del;
int_term = zeros(1, length(tspan));

%simulate system:
y = zeros(1, length(tspan));
y(1) = y_0;
int_term(1) = y_t - y_0;
for i = 1:(length(tspan) - 1)
    y(i + 1) = fdbck_PI_tank_lvl(y(i), y_t, d(i), k_c, tau_I, int_term(i), del);
    int_term(i + 1) = int_term(i) + del*(y_t - y(i + 1));
end

%plot tank level deviation vs time:
figure
plot(tspan, y)
xlabel("Time (min)", 'FontSize', 32)
ylabel("Y (tank ht deviation)", 'FontSize', 32)
title("PI Feedback Control with Measurement Error", 'FontSize', 36)
legend("Kc = " + k_c + ", Tau_I = " + tau_I/del + " min", 'FontSize', 26)
axis([0 600 -1 1])
ax = gca
ax.FontSize = 28;

function y_kp1 = fdbck_PI_tank_lvl(y_k, y_t, d_k, k_c, tau_I, I, del)
    %generate measured tank level deviation (with meas. error):
    y_m = y_k + 0.05*randn(1);
    %generate feedback control based on measured tank level deviation:
    u_k = k_c*(y_t - y_m) + (k_c/tau_I)*I;
    %compute tank level deviation at next time step:
    y_kp1 = y_k + del*(d_k - u_k);
end