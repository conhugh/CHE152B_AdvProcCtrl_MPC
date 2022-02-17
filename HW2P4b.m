%initialize parameters:
del = 1;   %sample time = 1 minute
tspan = [0:del:600];   %time span covers 600 min

%set initial condition & set point:
y_0 = 0;  
y_t = 0;

%initialize input vector:
d_1 = zeros(1, 300) + 0.5;
d_2 = zeros(1, 301) + 1;
d = cat(2, d_1, d_2);

%add measurement error:
dm = d + 0.05*randn(1, length(tspan));

%generate feedforward control input using control law  u(t) = dm(t):
u = dm;

%simulate system:
y = zeros(1, length(tspan));
y(1) = y_0;
for i = 1:(length(tspan) - 1)
    y(i + 1) = fdfwd_tank_lvl(y(i), d(i), u(i), del);
end

%plot tank level deviation vs time:
plot(tspan, y)
xlabel("Time (min)", 'FontSize', 32)
ylabel("Y (tank ht deviation)", 'FontSize', 32)
title("Feedforward Control with Measurement Error", 'FontSize', 36)
ax = gca;
ax.FontSize = 28;

function y_kp1 = fdfwd_tank_lvl(y_k, d_k, u_k, del)
    %compute tank level deviation at next time step:
    y_kp1 = y_k + del*(d_k - u_k);
end