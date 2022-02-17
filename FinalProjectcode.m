%Connor Hughes
%CH E 152B HW7


%% PID Tuning Revisited
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

K_c1 = 5*K_c1;
K_c2 = 6*K_c2;
Tau_I1 = 0.7*Tau_I1;
Tau_I2 = 1*Tau_I2;

tsam = 1;
dsys = c2d(ss(G), tsam);
ad = dsys.a;
bd = dsys.b;
cd = dsys.c;
%% Testing
close
G11 = tf(G.Numerator{1, 1}, G.Denominator{1, 1})
C1 = tf([K_c1, K_c1/Tau_I1], [1, 0])
% K_c1
% K_c1/Tau_I1
N11 = series(C1, G11)
CLsys11 = feedback(N11, 1)
zpk(CLsys11)

G12 = tf(G.Numerator{1, 2}, G.Denominator{1, 2})
C2 = tf([K_c2, K_c2/Tau_I2], [1, 0])
N12 = series(C2, G12)
CLsys12 = feedback(N12, 1)
zpk(CLsys12)
r = roots([1; 0.008729; 3.558e-05])
bode(CLsys12)

%% Simulation Test: 
tseries = linspace(1, 2000, 2000);
dnfreq = 0.0041;
T1_sp = zeros(1, length(tseries)) + 45;
u2_sp = zeros(1, length(tseries)) + 34.5;
for i = 1000:length(tseries)
    u2_sp(i) = 34.5 + 0.1*sin(tseries(i)/0.41);
end
T2_sp 
plot(tseries, T1_sp)
hold on
plot(tseries, u2_sp)
nsim = 2000
%% Simulation
%generate set points for closed-loop system:
%nsim = 3500;
% T1_sp = zeros(1, nsim) + 40;
% T1_sp(1001:1500) = 45;
% T1_sp(1501:2000) = 40;
% T1_sp(2001:2500) = 45;
% T1_sp(2501:end) = 40;
% T2_sp = zeros(1, nsim) + 40;
% T2_sp(1001:1500) = 50;
% T2_sp(1501:2000) = 45;
% T2_sp(2001:2500) = 40;
% T2_sp(2501:end) = 40;

T1_sp = T1_sp - Tstartavg(1);   %convert to deviation variables
u2_sp = u2_sp - Tstartavg(2);   %convert to deviation variables

%compute system response to above set points:
T_sp = [T1_sp; u2_sp]';
y0 = zeros(2, 1);
x0 = dsys.C\y0;
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
y_t = [T1_sp; u2_sp];
u(1, 1) = 30;
u(2, 1) = 30;

for k = 2:nsim
    y(:,k) = cd*x(:,k) + 0.0*randn(2,1);
    trackerr(:,k) = y_t(:,k)-y(:,k);
    u(1,k) = K_c1*(trackerr(1,k)+ 1/Tau_I1*interr(1,k)) + 30;
    %u(2,k) = K_c2*(trackerr(2,k)+ 1/Tau_I2*interr(2,k)) + 30;
    u(2, k) = y_t(2, k);
    if(k == nsim) 
        break 
    end
    interr(:,k+1) = interr(:,k)+trackerr(:,k)*tsam;
    x(:,k+1) = ad*x(:,k)+bd*(u(:,k) - 30);
end

y_t = y_t + Tstartavg';
y = y + Tstartavg';

time = linspace(1, nsim, nsim);
figure()
subplot(2,1,1)
ylabel("y", "rotation", 0)
plot(time, y(1, :), 'r', time, y(2, :), 'b', time, y_t, 'linewidth', 1.2)
ylim([35 55])
xlabel('time (s)')
ylabel('temperature (deg C)')
legend('Temperature 1', 'Temperature 2', 'Temp. 1 Set Point', 'Temp. 2 Set Point', 'FontSize', 12)
%title("P Gains: " + K_c1 + " " + K_c2 + ",  I Gains: " + Tau_I1 + " " + Tau_I2)
title("Temperatures and Set Points");
ax = gca
ax.FontSize = 16

subplot(2,1,2)
ylabel ("u", "rotation", 0)
stairs(time, u(1, :)', 'r', 'linewidth', 1.2);
hold on
stairs(time, u(2, :)', 'b', 'linewidth', 1.2);
xlabel('time (s)')
ylabel('Heater Setting (% Full Voltage)')
legend('Heater 1', 'Heater 2', 'FontSize', 12)
%title("P Gains: " + K_c1 + " " + K_c2 + ",  Taus: " + Tau_I1 + " " + Tau_I2)
title("Heater Settings")
ax = gca
ax.FontSize = 16

save('PIDvals', 'K_c1', 'K_c2', 'Tau_I1', 'Tau_I2')
%% Exercise 17: PID on TCLab Hardware
close all; clear all; clc
load PIDvals;
addpath('C:\Users\conno\Documents\UCSB Spring 2021 Homework\CH E 152B\tclab\MATLAB')

% include tclab.m for initialization
tclab;

disp('Test Heater 1')
disp('LED Indicates Temperature')

figure(1)
t1s = [];
t2s = [];
t1sp = [];
t2sp = [];
h1s = [];
h2s = [];

disp('Begin by running both heaters at 30%')
Tsp = [39.97; 36.93];
T1sp = 39.97;
T2sp = 36.93;
tsam = 1.0;
trackerr = [0; 0];
interr = [0; 0];
u1 = 30;
u2 = 30;

for i = 1:50
    tic;
    T1 = T1C();
    T2 = T2C();
    h1(u1);
    h2(u2);
    disp("Heater 1 setting: " + u1);
    disp("Heater 2 setting: " + u2);
    % LED brightness
    brightness1 = (T1 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness2 = (T2 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness = max(brightness1,brightness2);
    brightness = max(0,min(1,brightness)); % limit 0-1
    led(brightness);
    
    % plot heater and temperature data
    h1s = [h1s,u1];
    h2s = [h2s,u2];
    t1s = [t1s,T1];
    t2s = [t2s,T2];
    t1sp = [t1sp, T1sp];
    t2sp = [t2sp, T2sp];
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
    ylabel('Heater (% Duty Cycle)')
    xlabel('Time (sec)')
    legend('Heater 1','Heater 2','Location','NorthWest')
    drawnow;
    t = toc;
    interr = interr + trackerr*tsam;
    pause(max(0.01,1.0-t))
end
% for i = 1:500
%     tic;
%     T1 = T1C();
%     T2 = T2C();
%     trackerr= Tsp - [T1; T2];
%     u1 = K_c1*(trackerr(1) + 1/Tau_I1*interr(1));
%     u2 = K_c2*(trackerr(2) + 1/Tau_I2*interr(2));
%     h1(u1);
%     h2(u2);
%     disp("Heater 1 setting: " + u1)
%     disp("Heater 2 setting: " + u2)
%     % LED brightness
%     brightness1 = (T1 - 30)/50.0;  % <30degC off, >100degC full brightness
%     brightness2 = (T2 - 30)/50.0;  % <30degC off, >100degC full brightness
%     brightness = max(brightness1,brightness2);
%     brightness = max(0,min(1,brightness)); % limit 0-1
%     led(brightness);
%     
%     % plot heater and temperature data
%     h1s = [h1s,u1];
%     h2s = [h2s,u2];
%     t1s = [t1s,T1];
%     t2s = [t2s,T2];
%     n = length(t1s);
%     time = linspace(0,n+1,n);
%     clf
%     subplot(2,1,1)
%     plot(time,t1s,'r.','MarkerSize',10);
%     hold on
%     plot(time,t2s,'b.','MarkerSize',10);
%     ylabel('Temperature (degC)')
%     legend('Temperature 1','Temperature 2','Location','NorthWest')
%     subplot(2,1,2)
%     plot(time,h1s,'r-','LineWidth',2);
%     hold on
%     plot(time,h2s,'b--','LineWidth',2);
%     ylabel('Heater (% Duty Cycle)')
%     xlabel('Time (sec)')
%     legend('Heater 1','Heater 2','Location','NorthWest')
%     drawnow;
%     t = toc;
%     interr = interr + trackerr*tsam;
%     pause(max(0.01,1.0-t))
% end
for i = 1:3000
    tic;
    if i==5
        disp('Turn Temp 1 SP to 40, Temp 2 SP to 40')
        T1sp = 40;
        T2sp = 40;
    end
    if i==1005
        disp('Turn Temp 1 SP to 45, Temp 2 SP to 50')
        T1sp = 45;
        T2sp = 50;
    end
    if i==1505
        disp('Turn Temp 1 SP to 40, Temp 2 SP to 45')
        T1sp = 40;
        T2sp = 45;
    end
    if i==2005
        disp('Turn Temp 1 SP to 45, Temp 2 SP to 40')
        T1sp = 45;
        T2sp = 40;
    end
    if i==2505
        disp('Turn both set points to 40')
        T1sp = 40;
        T2sp = 40;
    end

    % read temperatures
    T1 = T1C();
    T2 = T2C();

    Tsp = [T1sp; T2sp];
    trackerr= Tsp - [T1; T2];
    u1 = K_c1*(trackerr(1) + 1/Tau_I1*interr(1));
    u2 = K_c2*(trackerr(2) + 1/Tau_I2*interr(2));
    disp('Heater 1 setting: ' + u1)
    disp('Heater 2 setting: ' + u2)
    h1(u1);
    h2(u2);
    
    % LED brightness
    brightness1 = (T1 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness2 = (T2 - 30)/50.0;  % <30degC off, >100degC full brightness
    brightness = max(brightness1,brightness2);
    brightness = max(0,min(1,brightness)); % limit 0-1
    led(brightness);
    
    % plot heater and temperature data
    h1s = [h1s,u1];
    h2s = [h2s,u2];
    t1s = [t1s,T1];
    t1sp = [t1sp, T1sp];
    t2sp = [t2sp, T2sp];
    t2s = [t2s,T2];
    n = length(t1s);
    time = linspace(0,n+1,n);
    clf
    subplot(2,1,1)
    plot(time,t1s,'r.','MarkerSize',10);
    hold on
    plot(time, t1sp, 'r--', 'MarkerSize',10);
    plot(time,t2s,'b.','MarkerSize',10);
    plot(time, t2sp, 'b--', 'MarkerSize',10);
    ylabel('Temperature (degC)')
    legend('Temperature 1', 'Temp 1 Set Point', 'Temperature 2', 'Temp 2 set point', 'Location','NorthWest')
    subplot(2,1,2)
    plot(time,h1s,'r-','LineWidth',2);
    hold on
    plot(time,h2s,'b--','LineWidth',2);
    ylabel('Heater (% Duty Cycle)')
    xlabel('Time (sec)')
    legend('Heater 1','Heater 2','Location','NorthWest')
    drawnow;
    t = toc;
    interr = interr + trackerr*tsam;
    pause(max(0.01,1.0-t))
end

disp('Turn off heaters')
h1(0);
h2(0);

disp('Heater Test Complete')

save('TCLabID_PID_1')

%% Re-Plot Results
figure
time = linspace(0, 2500, 2500);
subplot(2,1,1)
plot(time,t1s(551:end),'r.','MarkerSize',10);
hold on
plot(time,t2s(551:end),'b.','MarkerSize',10);

plot(time, t1sp(551:end), time, t2sp(551:end), 'linewidth', 1.2)
ylabel('Temperature (degC)')
legend('Temperature 1','Temperature 2', 'Temp 1 Set Point', 'Temp 2 Set Point', 'Location','NorthEast', 'FontSize', 12)
ax = gca
ax.FontSize = 16;

subplot(2,1,2)
plot(time,h1s(551:end),'r-','LineWidth',2);
hold on
plot(time,h2s(551:end),'b--','LineWidth',2);
ylabel('Heater (% Duty Cycle)')
xlabel('Time (sec)')
legend('Heater 1','Heater 2','Location','NorthWest', 'FontSize', 12)
ax = gca
ax.FontSize = 16;

    
    