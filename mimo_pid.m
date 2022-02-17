% 
%  MIMO PID with various pairings
% 
%  jbr, 12/2/99
%  updated 3/31/11
%  updated 2/17/12
%

clear all
close all

data = struct();

k11   = 2; tau11 = 10;
k12   = 2; tau12 = 1;
k21   = 1; tau21 = 1;
k22   = -4; tau22 = 10;
K=[k11 k12; k21 k22];
RGA = K .* inv(K)';

tsam = 0.1;
num = {k11, k12; k21, k22};
den = {[tau11 1], [tau12 1]; [tau21 1], [tau22 1]};
sys = tf(num, den);

dsys = c2d(ss(tf(num,den)), tsam);
ad = dsys.a;
bd = dsys.b;
cd = dsys.c;

pid_cases = {'rga', 'dyn'};
for ncase = 1: length(pid_cases)
  part = pid_cases{ncase};
  if (strcmp(part, 'rga'))
    %% RGA pairing
    k_c1   = 2;
    tau_i1 = 2;
    k_c2   = -2;
    tau_i2 = 3;
    u1y = 1;
    u2y = 2;
  elseif (strcmp(part, 'dyn'))
    %% dynamic pairing; reverse pairing
    k_c1   = 5;
    tau_i1 = 2;
    k_c2   = 5;
    tau_i2 = 3;
    u1y = 2;
    u2y = 1;
  else
    error ('PID case error')
  end%if
  %%
  %% pid simulation
  %%
  nsim=max([tau11,tau12,tau21,tau22])/tsam;
  time=(0:nsim-1)*tsam;
  order = length(ad(:,1));
  x = zeros(order,nsim);
  interr = zeros(2,nsim);
  y = zeros(2, nsim);
  y_t = [ones(1,nsim); zeros(1,nsim)];
  y_t(:,1) = 0;

  for k = 2:nsim
    y(:,k) = cd*x(:,k) + 0.0*randn(2,1);
    trackerr(:,k) = y_t(:,k)-y(:,k);
    u(1,k) = k_c1*( trackerr(u1y,k)+ 1/tau_i1*interr(u1y,k) );
    u(2,k) = k_c2*( trackerr(u2y,k)+ 1/tau_i2*interr(u2y,k) );
    if(k == nsim) break endif
    interr(:,k+1) = interr(:,k)+trackerr(:,k)*tsam;
    x(:,k+1) = ad*x(:,k)+bd*u(:,k);
  endfor
  data.(part) = [time' u' y' y_t'];
  figure()
  subplot(2,1,1)
  ylabel("y", "rotation", 0)
  plot(time, y, time, y_t)
  subplot(2,1,2)
  ylabel ("u", "rotation", 0)
  stairs(time, u')
endfor

gnuplotsave('mimo_pid.dat', data)