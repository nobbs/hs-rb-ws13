% RB-System aufstellen
Q_a   = size(Ak, 1);
Ak_rb = cell(Q_a, 1);
for i = 1:Q_a
	Ak_rb{i} = Z' * Ak{i} * Z;
end
f_rb = Z' * f;

% Vorbereitung für Fehlerberechnung
G = prepare_error_estimator(X, Ak, f, Z);

% Zufällig Parameter wählen
random_mu = @() generate_parameter_grid(mu_min, mu_max, P, 1, 'mc_log');
