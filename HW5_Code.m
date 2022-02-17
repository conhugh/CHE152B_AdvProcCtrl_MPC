%Connor Hughes
%CH E 152B
%HW 5

%% Exercise 11
clear, clc, close all;
load data-exercise-siso-1
figure
subplot(2, 1, 1)
plot(tsim, u, 'g')
ax = gca;
ax.FontSize = 20;
xlabel('time')
ylabel('u')
ylim([-1.2, 1.2])
lgd = legend('Input Signal');
lgd.FontSize = 14;
subplot(2, 1, 2)
scatter(tsim, ym, 'Marker', '.', 'MarkerEdgeColor', 'r')
ax = gca;
ax.FontSize = 20;
xlabel('time')
ylabel('y')
hold on

dat = iddata(ym', u', 0.1);
% sys = tfest(dat, 1);
%[pv1, pvsd1] = getpvec(sys)
sys = tfest(dat, 2);
[pv2, pvsd2] = getpvec(sys)
% sys = tfest(dat, 3);
% [pv3, pvsd3] = getpvec(sys)
% sys = tfest(dat, 4);
% [pv4, pvsd4] = getpvec(sys)

yest = lsim(sys, u, tsim)
plot(tsim, yest, 'b')
lgd = legend('Measured Response', 'Model Response');
lgd.FontSize = 14;

figure
sys_tf = tf(sys.Numerator, sys.Denominator);
step(sys_tf)
ax = gca;
ax.FontSize = 20;
ylabel('Amplitude', 'FontSize', 20)
xlabel('Time', 'FontSize', 20)
title('Step Response', 'FontSize', 20)

%% Exercise 12
clear, clc, close all;
load data-exercise-siso-2
figure
subplot(2, 1, 1)
plot(tsim, u, 'g')
ax = gca;
ax.FontSize = 20;
xlabel('time')
ylabel('u')
ylim([-1.2, 1.2])
lgd = legend('Input Signal');
lgd.FontSize = 14;
subplot(2, 1, 2)
scatter(tsim, ym, 'Marker', '.', 'MarkerEdgeColor', 'r')
ax = gca;
ax.FontSize = 20;
xlabel('time')
ylabel('y')
hold on

dat = iddata(ym', u', 0.5);
sys = tfest(dat, 1, 0, NaN);
[pv1, pvsd1] = getpvec(sys)
% sys = tfest(dat, 2);
% [pv2, pvsd2] = getpvec(sys)
% sys = tfest(dat, 3);
% [pv3, pvsd3] = getpvec(sys)
% sys = tfest(dat, 4);
% [pv4, pvsd4] = getpvec(sys)

yest = lsim(sys, u, tsim);
plot(tsim, yest, 'b')
lgd = legend('Measured Response', 'Model Response');
lgd.FontSize = 14;

figure
sys_tf = tf(sys.Numerator, sys.Denominator);
step(sys_tf)
ax = gca;
ax.FontSize = 20;
ylabel('Amplitude', 'FontSize', 20)
xlabel('Time', 'FontSize', 20)
title('Step Response', 'FontSize', 20)

%% Exercise 13
clear, clc, close all;
load data-exercise-mimo.mat
figure
subplot(2, 2, 1)
plot(tsim, u(1, :), 'g')
ylabel('u_1')
xlabel('time')
ylim([-1.2, 1.2])
ax = gca
ax.FontSize = 20;
subplot(2, 2, 3)
plot(tsim, u(2, :), 'g')
ylabel('u_2')
xlabel('time')
ylim([-1.2, 1.2])
ax = gca
ax.FontSize = 20;
subplot(2, 2, 2)
scatter(tsim, ym(1, :), 'Marker', '.', 'MarkerEdgeColor', 'r')
ylabel('y_1')
xlabel('time')
hold on
ax = gca
ax.FontSize = 20;
subplot(2, 2, 4)
scatter(tsim, ym(2, :), 'Marker', '.', 'MarkerEdgeColor', 'r')
ylabel('y_2')
xlabel('time')
ax = gca
ax.FontSize = 20;
hold on

id = iddata(ym', u', 0.1)
Gsys = tfest(id, 1)


%figure
G = ssest(id, 2);
L = lsim(Gsys, u, tsim);
subplot(2, 2, 2)
plot(tsim, L(:, 1), 'b')
subplot(2, 2, 4)
plot(tsim, L(:, 2), 'b')

figure
step(Gsys)



