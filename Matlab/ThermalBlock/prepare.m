%% Mesh erzeugen und verfeinern
[p, e, t] = initmesh(geom, 'Hmax', 0.075);
% [p, e, t] = initmesh(geom, 'Hmax', 0.05);
% [p, e, t] = refinemesh(geom, p, e, t);

%% Parameterbereich festlegen
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
N_max = 100;
tolerance = 1e-8;

% Trainings-Parameter
Xi_train = generate_parameter_grid(mu_min, mu_max, P, n_train, 'mc_lin');
n_train = size(Xi_train, 2);

% Test-Parameter
if tgl_train_eq_test == 1
    Xi_test = Xi_train;
else
    Xi_test = generate_parameter_grid(mu_min, mu_max, P, n_test, 'mc_log');
end
n_test = size(Xi_test, 2);
