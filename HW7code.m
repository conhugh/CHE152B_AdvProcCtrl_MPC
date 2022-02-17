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
nsim = 6000;
T1_sp = zeros(1, nsim) + 40;
T1_sp(1501:2500) = 45;
T1_sp(2501:3500) = 40;
T1_sp(3501:4500) = 45;
T1_sp(4501:end) = 40;
T2_sp = zeros(1, nsim) + 40;
T2_sp(1501:2500) = 50;
T2_sp(2501:3500) = 45;
T2_sp(3501:4500) = 40;
T2_sp(4501:end) = 40;

T1_sp = T1_sp - Tstartavg(1);   %convert to deviation variables
T2_sp = T2_sp - Tstartavg(2);   %convert to deviation variables

%compute system response to above set points:
T_sp = [T1_sp; T2_sp]';
y0 = zeros(2, 1);
x0 = dsys.C\y0;
% time=(0:nsim-1)*tsam;
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
u(1, 1) = 30;
u(2, 1) = 30;

for k = 2:nsim
    y(:,k) = cd*x(:,k) + 0.0*randn(2,1);
    trackerr(:,k) = y_t(:,k)-y(:,k);
    u(1,k) = K_c1*(trackerr(1,k)+ 1/Tau_I1*interr(1,k)) + 30;
    u(2,k) = K_c2*(trackerr(2,k)+ 1/Tau_I2*interr(2,k)) + 30;
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

%% Exercise 17: PID on TCLab Hardware
close all; clear all; clc
load PIDvals;

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

disp('Begin by driving both temperature readings to 40 deg C')
T1sp = 40;
T2sp = 40;
Tsp = [40; 40];
tsam = 1.0;
trackerr = [0; 0];
interr = [0; 0];


for i = 1:1000
    tic;
    T1 = T1C();
    T2 = T2C();
    trackerr= Tsp - [T1; T2];
    u1 = K_c1*(trackerr(1) + 1/Tau_I1*interr(1));
    u2 = K_c2*(trackerr(2) + 1/Tau_I2*interr(2));
    h1(u1);
    h2(u2);
    disp("Heater 1 setting: " + u1)
    disp("Heater 2 setting: " + u2)
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
for i = 1:3000
    tic;
    if i==505
        disp('Turn Temp 1 SP to 85, Temp 2 SP to 35')
        T1sp = 85;
        T2sp = 35;
    end
%     if i==1505
%         disp('Turn Temp 1 SP to 40, Temp 2 SP to 45')
%         T1sp = 40;
%         T2sp = 45;
%     end
%     if i==2505
%         disp('Turn Temp 1 SP to 45, Temp 2 SP to 40')
%         T1sp = 45;
%         T2sp = 40;
%     end
%     if i==3505
%         disp('Turn both set points to 40')
%         T1sp = 40;
%         T2sp = 40;
%     end

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
    t2s = [t2s,T2];
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

disp('Turn off heaters')
h1(0);
h2(0);

disp('Heater Test Complete')

save('TCLabID4')

%% Re-Plot Results
figure
time = linspace(1001, 6000, 5000);
subplot(2,1,1)
plot(time,t1s(1001:end),'r.','MarkerSize',10);
xlim([1000, 6000])
hold on
plot(time,t2s(1001:end),'b.','MarkerSize',10);
T1spHist = zeros(1, 5000);
T2spHist = zeros(1, 5000);
T1spHist(1:500) = T1spHist(1:500) + 40;
T2spHist(1:500) = T2spHist(1:500) + 40;
T1spHist(501:1500) = T1spHist(501:1500) + 45;
T2spHist(501:1500) = T2spHist(501:1500) + 50;
T1spHist(1501:2500) = T1spHist(1501:2500) + 40;
T2spHist(1501:2500) = T2spHist(1501:2500) + 45;
T1spHist(2501:3500) = T1spHist(2501:3500) + 45;
T2spHist(2501:3500) = T2spHist(2501:3500) + 40;
T1spHist(3501:5000) = T1spHist(3501:5000) + 40;
T2spHist(3501:5000) = T2spHist(3501:5000) + 40;
plot(time, T1spHist, time, T2spHist, 'linewidth', 1.2)
ylabel('Temperature (degC)')
legend('Temperature 1','Temperature 2', 'Temp 1 Set Point', 'Temp 2 Set Point', 'Location','NorthEast', 'FontSize', 12)
ax = gca
ax.FontSize = 16;

subplot(2,1,2)
plot(time,h1s(1001:end),'r-','LineWidth',2);
hold on
plot(time,h2s(1001:end),'b--','LineWidth',2);
xlim([1000, 6000])
ylabel('Heater (% Duty Cycle)')
xlabel('Time (sec)')
legend('Heater 1','Heater 2','Location','NorthWest', 'FontSize', 12)
ax = gca
ax.FontSize = 16;

    
    