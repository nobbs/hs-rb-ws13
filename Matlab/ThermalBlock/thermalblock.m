%% 3x3 Thermal-Domain laden
load('dreixdrei_geom.mat');

%% Mesh erzeugen und verfeinern
[p, e, t] = initmesh(geom);
[p, e, t] = refinemesh(geom, p, e, t);
% [p, e, t] = refinemesh(geom, p, e, t);
% [p, e, t] = refinemesh(geom, p, e, t);

%% Parameterbereich festlegen
bx = 3;
by = 3;
P = bx * by - 1;

mu_r = 10;
mu_min = 1 / sqrt(mu_r);
mu_max = sqrt(mu_r);
mu_log_min = log(mu_min);
mu_log_max = log(mu_max);
Mu_bar = ones(P, 1);

% Parameter-Bereiche
D_box = repmat([mu_min, mu_max], P, 1);
D_ln_box = repmat([mu_log_min, mu_log_max], P, 1);

%% Offline-Stage
n_train = 1000;
Xi = generate_parameter_grid(mu_min, mu_max, P, n_train, 'mc_lin');
Z = reduced_basis(p, e, t, bx, by, Xi, size(p, 2), Mu_bar);

%% Online-Stage
Mu_test = generate_parameter_grid(mu_min, mu_max, P, 1, 'mc_lin');
[A, Ak, F, X] = assemble(p, e, t, bx, by, Mu_test);
U_rb = rb_online(Mu_test, Z, Ak, F);
S_rb = F' * U_rb

U_fe = A \ F;
S_fe = F' * U_fe

pdesurf(p, t, abs(U_rb - U_fe));
xlabel('x');
ylabel('y');
colormap jet;

% [A, Ak, f, X] = assemble(p, e, t, bx, by, Mu_bar);
% U = A \ f;

% Plotten
% pdesurf(p, t, U);
% colormap jet;
% view(0, 90)
