% alpha_LB = @(mu) min([mu / mu_bar, 1 / mu_bar])
% gamma_UB = @(mu) max([mu / mu_bar, 1 / mu_bar])

% Provisorische Fehlerüberprüfung
% Errors = zeros(length(G_lin_aq), 1);
% for i = 1:length(Errors)
% 	mu = G_lin_aq(i);
% 	A_rb = Z' * (mu * Ak(:, :, 1) + Ak(:, :, 2)) * Z;
% 	F_rb = Z' * F;
% 	U_rb = Z * (A_rb \ F_rb);
% 	U_ex = exakteLoesung(Grid, mu);
% 	Errors(i) = max(abs(U_rb - U_ex));
% end

% plot(G_lin_aq, Errors);

% Nach Quelle:

H = zeros(n, 1 + 2 * N_max);
H(:, 1) = F;
H(:, 2:2:end) = -Ak(:, :, 1)' * Z;
H(:, 3:2:end) = -Ak(:, :, 2)' * Z;

G = H' * (X \ H);
