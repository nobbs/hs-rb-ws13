A_rb = Z' * (mu_online * Ak(:, :, 1) + Ak(:, :, 2)) * Z;
F_rb = Z' * F;
U_rb = Z * (A_rb \ F_rb);
U_ex = exakteLoesung(Grid, mu_online);
plot(Grid, U_rb - U_ex);
