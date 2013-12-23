% Lösen
A_rb_check = Z' * (mu_check * Ak(:, :, 1) + Ak(:, :, 2)) * Z;
F_rb_check = Z' * F;
U_rb_check = A_rb_check \ F_rb_check;

% Fehler berechnen
Eps = zeros(1 + 2 * N_max, 1);
Eps(1) = 1;
Eps(2:2:end) = mu_check * U_rb_check;
Eps(3:2:end) = U_rb_check;

% Fehler
resid = sqrt(Eps' * G * Eps)
alpha_LB = min([mu_check / mu_bar, 1 / mu_bar])
Err_en = resid / sqrt(alpha_LB)

