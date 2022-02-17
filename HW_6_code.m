%Connor Hughes
%CH E 152B HW 6
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

%allow time to warm up and reach steady state
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
%begin pulse test:
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
U = U - 30;   %SWITCH TO DEVIATION VARIABLES
Tendavg = mean(T(4700:end, :));
Tstartavg = mean(T(200:600, :));
tr1 = linspace(Tstartavg(1), Tendavg(1), 5000);
tr2 = linspace(Tstartavg(2), Tendavg(2), 5000);
tr1 = tr1 - mean(tr1);
tr2 = tr2 - mean(tr2);

Tdt(:, 1) = T(:, 1) - tr1';
Tdt(:, 2) = T(:, 2) - tr2';

Tdt = Tdt - Tstartavg; %SWITCH TO DEVIATION VARIABLES

%plot data (with small linear trend removed)
figure
subplot(2, 1, 1)
plot(time,Tdt(:, 1),'r.','MarkerSize',10);
hold on
plot(time,Tdt(:, 2),'b.','MarkerSize',10);
ylabel('Temperature (degC)')
subplot(2,1,2)
plot(time,h1s - 30,'r-','LineWidth',2);   %PLOT DEV VARS
hold on
plot(time,h2s - 30,'b--','LineWidth',2);  %PLOT DEV VARS
ylabel('Heater (% duty cycle - 30)')
xlabel('Time (sec)')
legend('Heater 1','Heater 2','Location','NorthWest')
ax = gca
ax.FontSize = 18
drawnow;

% %use system ID toolbox functions to estimate system from data
% dat = iddata(Tdt(600:end, :), U(600:end, :), 1.0004);
% tsim = time(600:end);
% [G, ic] = tfest(dat, 1);
% [pv, pvsd] = getpvec(G)
% Gss = ssest(dat, 5);
% %Gssd = c2d(Gss, 1.0004);
% Gssd = c2d(ss(G), 1.0004)

% %simulate solution of estimated system for the same input used in the pulse
% %test and plot on the same axes as the test data, for comparison
% subplot(2, 1, 1)
% x0 = Gssd.C\(Tdt(600, :)');
% y = lsim(Gssd, U(600:end, :), tsim, x0, 'zoh');
% plot(tsim, y, 'g', 'LineWidth', 1)
% legend('Temp. 1 Detrended','Temp. 2 Detrended', 'Temp. 1 Simulated', 'Temp. 2 Simulated', 'Location','NorthWest')
% save('TCLabID4comp')

dat = iddata(Tdt, U, 1.0004);
tsim = time;
[G, ic] = tfest(dat, 1);
[pv, pvsd] = getpvec(G)
Gss = ssest(dat, 5);
%Gssd = c2d(Gss, 1.0004);
Gssd = c2d(ss(G), 1.0004)

subplot(2, 1, 1)
x0 = Gssd.C\(Tdt(1, :)');
y = lsim(Gssd, U, tsim, x0, 'zoh');
plot(tsim, y, 'g', 'LineWidth', 1)
ax = gca
ax.FontSize = 18
legend('Temp. 1 Detrended','Temp. 2 Detrended', 'Temp. 1 Simulated', 'Temp. 2 Simulated', 'Location','NorthWest')
title('Pulse Test Data in deviation variables')
save('TCLabID4comp')


%% TF conversion
close all
Gssd = c2d(ss(G), 1.004)
figure
step(G)

%% PID Tuning
close all;
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
y0 = [28; 28];
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
u(1, 1) = 0;
u(2, 1) = 0;

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

time = linspace(1, 4500, 4500);
figure()
subplot(2,1,1)
ylabel("y", "rotation", 0)
plot(time, y(1, 13501:end), 'r', time, y(2, 13501:end), 'b', time, y_t(:, 13501:end), 'linewidth', 1.2)
ylim([30 60])
xlabel('time (s)')
ylabel('temperature (deg C)')
legend('Temperature 1', 'Temperature 2', 'Temp. 1 Set Point', 'Temp. 2 Set Point', 'FontSize', 12)
title("Temperatures and Set Points");
ax = gca
ax.FontSize = 16

subplot(2,1,2)
ylabel ("u", "rotation", 0)
stairs(time, u(1, 13501:end)', 'r', 'linewidth', 1.2);
hold on
stairs(time, u(2, 13501:end)', 'b', 'linewidth', 1.2);
xlabel('time (s)')
ylabel('Heater Setting (% Full Voltage)')
legend('Heater 1', 'Heater 2', 'FontSize', 12)
title("Heater Settings")
ax = gca
ax.FontSize = 16

%% Step tests 
step(Gssd)
step(G)

%% Plot: 
plot(time, t1s, time, t2s)

