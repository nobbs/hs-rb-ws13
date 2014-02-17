%% Mesh erzeugen und ggf. verfeinern
[p, e, t] = initmesh(geom, 'Hmax', 0.075);
% [p, e, t] = refinemesh(geom, p, e, t);

%% Parameterbereich
% Anzahl Parameter
P = bx * by - 1;

% Parameterbereich festlegen. (P > 1 => Würfel)
mu_r = 100;
mu_min = 1 / sqrt(mu_r);
mu_max = sqrt(mu_r);

% Referenzparameter
mu_bar = ones(P, 1);

%% FE-Steifigkeitsmatrizen assemblieren
[Ak, f, B, Ud] = assemble_fe(p, e, t, bx, by);
% Hilfsfunktion, um Dirichlet-RB. wieder reinzubringen
assemble_sol = @(u) B * u + Ud;

% Matrix der diskreten X-Norm
X = assemble_fe_A(Ak, [mu_bar; 1]);

%% Offline-Stage Daten
% Maximale RB-Größe
N_max = 100;
% Zu erreichende Fehler-Toleranz
tolerance = 1e-8;

% Trainings-Parameter erzeugen
Xi_train = generate_parameter_grid(mu_min, mu_max, P, n_train, 'mc_log');
n_train = size(Xi_train, 2);

% Test-Parameter erzeugen (oder Trainings-Parameter benutzen)
if tgl_train_eq_test == 1
    Xi_test = Xi_train;
else
    Xi_test = generate_parameter_grid(mu_min, mu_max, P, n_test, 'mc_log');
end
n_test = size(Xi_test, 2);
