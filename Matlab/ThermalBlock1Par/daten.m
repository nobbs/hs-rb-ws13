%% 2x1 Thermal-Domain laden
load('geom_zweixeins.mat');

%% Mesh erzeugen und verfeinern
[p, e, t] = initmesh(geom, 'Hmax', 0.085);
% [p, e, t] = refinemesh(geom, p, e, t);

%% Parameterbereich festlegen
bx = 2;
by = 1;
P = bx * by - 1;

mu_r = 100;
mu_min = 1 / sqrt(mu_r);
mu_max = sqrt(mu_r);
mu_log_min = log(mu_min);
mu_log_max = log(mu_max);
mu_bar = ones(P, 1);

% Parameter-Bereiche
% D_box = repmat([mu_min, mu_max], P, 1);
% D_ln_box = repmat([mu_log_min, mu_log_max], P, 1);

%% Steifigkeitsmatrizen assemblieren
[Ak, F, B, Ud] = fe_assemble(p, e, t, bx, by);
% Hilfsfunktion, um Dirichlet-RB. wieder reinzubringen
assemble_sol = @(u) B * u + Ud;

% Matrix der diskreten X-Norm
X = fe_assemble_A(Ak, mu_bar);

num_ps_grid = size(p, 2);
num_ps_no_dbc = size(Ak{1}, 1);

%% Offline-Stage Daten
N_max = 10;
tolerance = 1e-6;

% Trainings-Parameter
n_train = 1e4;
Xi_train = generate_parameter_grid(mu_min, mu_max, P, n_train, 'mc_log');
% Test-Parameter
n_test = 1e4;
Xi_test = generate_parameter_grid(mu_min, mu_max, P, n_test, 'mc_log');
