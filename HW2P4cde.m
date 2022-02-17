%% Problem 4c,d, and e: Feedback Control (proportional only)
%initialize proportional gain:
k_c = -2;

%simulate system:
y = zeros(1, length(tspan));
y(1) = y_0;
for i = 1:(length(tspan) - 1)
    y(i + 1) = fdbck_P_tank_lvl(y(i), y_t, d(i), k_c, del);
end

%plot tank level deviation vs time:
plot(tspan, y)
xlabel("Time (min)", 'FontSize', 32)
ylabel("Y (tank ht deviation)", 'FontSize', 32)
title("Proportional Feedback Control with Measurement Error", 'FontSize', 36)
legend("Kc = -2", 'FontSize', 26)
%legend("Kc = -1")
%legend("Kc = -2")
ax = gca;
ax.FontSize = 28;

function y_kp1 = fdbck_P_tank_lvl(y_k, y_t, d_k, k_c, del)
    %generate measured tank level deviation (with meas. error):
    y_m = y_k + 0.05*randn(1);
    %generate feedback control based on measured tank level deviation:
    u_k = k_c*(y_t - y_m);
    %compute tank level deviation at next time step:
    y_kp1 = y_k + del*(d_k - u_k);
end