% RB-Gls aufstellen
num_Ak = size(Ak, 1);
Ak_rb = cell(num_Ak, 1);
for i = 1:num_Ak
	Ak_rb{i} = Z' * Ak{i} * Z;
end
F_rb = Z' * F;

% Vorbereitung für Fehlerberechnung
if tgl_certify == 1
    estimate_error_prep;
end

% Zufällig Parameter wählen
random_mu = @() generate_parameter_grid(mu_min, mu_max, P, 1, 'mc_lin');
