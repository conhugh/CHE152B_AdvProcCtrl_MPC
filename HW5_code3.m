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