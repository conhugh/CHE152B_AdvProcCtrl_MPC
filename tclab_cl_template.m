import casadi.*
mpc = import_mpctools();

figure(1)
t1s = [];
t2s = [];
h1s = [];
h2s = [];
yt1s = [];
yt2s = [];
% initial heater values
ht1 = xx;
ht2 = xx;

%% mode of operation;
%% either testing in simulation with nominal model
%% or applying to the tclab/arduino plant hardware

mode = 'simulation'
%%mode = 'plant'

%% use either R or S for tuning the MPC regulator
tune = 'R'
%%tune = 'S'

%% use either input or output integrating disturbance model to remove offset
%%dmodel = 'input'
dmodel = 'output'

if (strcmp(mode, 'plant'))
  h1(ht1);
  h2(ht2);
end%if


%% insert identified model here

A = [];
B = [];
C = [];
D = [];

[p, n] = size(C);
[n, m] = size(B);

%% choose the steady state us0 used in ID experiment,
us0 = [ht1, ht2]';
%% corresponding temperature steady states, ys
ys0 = [xx, xx]';

%% start plant simulation state at steady state; xs0=0  
xs0 = zeros(n,1);

K = C*inv(eye(n)-A)*B;
RGA = K.*inv(K)';

Qy = diag([1,1])*diag([1./(ys0.*ys0)]);
Q = C'*Qy*C;

switch tune
  case 'R'
    R = 0.6*diag([1./(us0.*us0)]);
    S = zeros(m);
  case 'S'
    R = zeros(m);
    S = 5*diag([1./(us0.*us0)]);
  otherwise
    error('Unknown choice for tuning: %s', tune);
end%switch

switch dmodel
  case 'output'
    Bd = zeros(n,m);
    Cd = eye(p);
    Nd = p;
  case 'input'
    Bd = B;
    Cd = zeros(p, m);
    Nd = m;
  otherwise
    error('Unknown choice for disturbance model: %s', dmodel);
end%switch
Hor = 50;
Nx = n;
Nu = m;
Aaug = [A, Bd; zeros(Nd, Nx), eye(Nd)];
Baug = [B; zeros(Nd, Nu)];
Caug = [C, Cd];
Naug = size(Aaug,1);
%% Detectability test of disturbance model
detec = rank([eye(Naug) - Aaug; Caug]);
if detec < Nx + Nd
  fprintf(' * Augmented system is not detectable!\n')
end%if
small = 1e-5; % Small number.
Qw = zeros(Naug);
Qw(1:Nx, 1:Nx) = small*eye(Nx);
Qw(Nx+1:end, Nx+1:end) = eye(Nd);
Rv = eye(p);
[L, ~, Pe] = dlqe(Aaug, eye(Naug), Caug, Qw, Rv);
Lx = L(1:Nx,:);
Ld = L(Nx+1:end,:);

%% heater constraints
lb = struct();
ub = struct();
lb.u = zeros(m,1);
ub.u = 100*ones(m,1);
lb.Du = -inf*(ub.u-lb.u);
ub.Du = - lb.Du;

% ## build up casadi mpc functions 
f = @(x, u, xs0, us0) xs0 + A*(x-xs0) + B*(u-us0);
N = struct('x', n, 'u', m, 't', Hor);

% ## create MPC regulator
fcasadi = mpc.getCasadiFunc(f, [N.x, N.u, N.x, N.u], ...
			    {'x', 'u', 'xs0', 'us0'}, {'f'});
parreg = struct('xs', zeros(n,1), 'us', us0, 'uprev', us0, ...
		'xs0', xs0, 'us0', us0);
l = @(x, u, Du, xs, us) (x-xs)'*Q*(x-xs) + (u-us)'*R*(u-us) + Du'*S*Du;
lcasadi = mpc.getCasadiFunc(l, [N.x, N.u, N.u, N.x, N.u], {'x', 'u', 'Du', 'xs', 'us'}, {'l'});
regulator = mpc.nmpc('f', fcasadi, 'N', N, 'x0', zeros(n,1), 'l', lcasadi, ...
		     'verbosity', 1, 'lb', lb, 'ub', ub, ...
		     'par', parreg);
% ## create ss target optimizer
fss = @(x, u, ds, xs0, us0) xs0 + A*(x-xs0) + B*(u-us0) +Bd*ds;
fsscasadi = mpc.getCasadiFunc(fss, [N.x, N.u, Nd, N.x, N.u], ...
			      {'x', 'u', 'ds', 'xs0', 'us0'}, {'fss'});
lss = @(x, u, ds, ysp, ys0) (ysp - (C*x + Cd*ds + ys0))'*Qy *(ysp - (C*x + Cd*ds + ys0)); 
lsscasadi = mpc.getCasadiFunc(lss, [N.x, N.u, Nd, p, p], {'x', 'u', 'ds', 'ysp', 'ys0'}, {'lss'});
parss = struct('ds', zeros(p,1), 'ysp', zeros(p,1), 'xs0', xs0, 'us0', us0, 'ys0', ys0);
Nss = struct('x', n, 'u', m);
sstarg = mpc.sstarg('f', fsscasadi, 'N', Nss, 'l', lsscasadi, 'lb', lb, ...
		    'ub', ub, 'par', parss); 

%% setpoint trajectory
ytseq = [ys0, ys0 + [5; -5], ys0, ys0 + [-5; 5]];
%%ytseq = [ys0, ys0, ys0, ys0];
distseq  = 0*[[1;0], [0;0], [0; -1], [0;0]];


%% mpc simulation
nmin = 5;
[~, nsp] = size(ytseq);
nts = nsp*nmin*60;
yp = NaN(p, nts);
xp = NaN(n, nts);
xp(:, 1) = xs0;
xs = xp;
xhat = xp;
xhatm = zeros(n,1);



dhat = NaN(p, nts);
dhatm = zeros(p,1);
us = NaN(m, nts);

yvar = 0.25;
it = 0;

for j = 1:nsp
  yt = ytseq(:,j);
  dist = distseq(:, j);
  yt1 = yt(1);
  yt2 = yt(2);
%% run to steady state with this setpoint
for i = 1:nmin*60
  tic;
  it = it + 1;
  if (strcmp(mode, 'plant'))
    %% read temperatures
    t1 = T1C();
    t2 = T2C();
  else
    %% simulate temperatures
    ytmp = C*xp(:, it) + ys0 + sqrt(yvar)*randn(p,1) + Cd*dist;
    t1 = ytmp(1);
    t2 = ytmp(2);
  end%if
  yp(:, it) = [t1;t2];
  ey = yp(:, it) - (C*xhatm + Cd*dhatm + ys0);
  xhat(:, it) = xhatm + Lx*ey;
  dhat(:, it) = dhatm + Ld*ey;
  %% compute steady-state targets
  sstarg.par.ysp = yt;
  sstarg.par.ds = dhat(:, it);
  sstarg.solve();
  xs(:, it) = sstarg.var.x;
  us(:, it) = sstarg.var.u;

  %% compute optimal dynamic control
  x0 = xhat(:, it);
  regulator.fixvar('x', 1, x0);
  regulator.par.xs = xs(:, it);
  regulator.par.us = us(:, it);
  regulator.solve();
  u(:, it) = regulator.var.u(:,1);
  ht1 = u(1, it);
  ht2 = u(2, it);
  if (strcmp(mode, 'plant'))
    %% send heater settings to board
    h1(ht1);
    h2(ht2);
  else
    %% test in simulation mode; send heater settings to model
    xp(:, it+1) = f(xp(:, it), u(:, it), xs0, us0) + Bd*dist;
  end%if
  %% update previous input and advance state estimate
  regulator.par.uprev = u(:, it);
  %% advance state of controller model and estimated disturbance
  xhatm = f(xhat(:, it), u(:, it), xs0, us0) + Bd*dhat(:, it); 
  dhatm = dhat(:, it);
  %% plot heater and temperature data
  h1s = [h1s, ht1];
  h2s = [h2s, ht2];
  t1s = [t1s, t1];
  t2s = [t2s, t2];
  yt1s = [yt1s, yt1];
  yt2s = [yt2s, yt2];
  nplot = length(t1s);
  time = linspace(0, nplot+1, nplot);
  clf
  subplot(2,1,1)
  hold on
  plot(time,t1s,'r.','MarkerSize',10)
  plot(time, yt1s, 'r--');
  plot(time, t2s,'b.','MarkerSize',10)
  plot(time, yt2s, 'b--');
  ylabel('Temperature (degC)')
  legend('Temperature 1', 'T1sp', 'Temperature 2', 'T2sp', 'Location', 'NorthWest')
  drawnow;
  subplot(2,1,2)
  hold on
  stairs(time,h1s,'r-','LineWidth',2);
  stairs(time,h2s,'b--','LineWidth',2);
  ylabel('Heater (0-5.5 V)')
  xlabel('Time (sec)')
  legend('Heater 1','Heater 2','Location','NorthWest')
  drawnow;
  t = toc;
  if (t >= 1)
    disp ('warning, sample time t = ')
    disp (t)
  end%if
  if (strcmp(mode, 'plant'))
    %% hold for sample time of 1 sec
    pause(max(0.01,1.0-t))
  end%if
end%for
%disp ('update heater power')
%keyboard
end%for

table = [time; h1s; h2s; t1s; t2s; yt1s; yt2s; dhat]';
save 'tclab_cl_template.dat' table