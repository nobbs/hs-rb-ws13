% Theta aus dem gewählten Parameter konstruieren
Theta = [online_mu; 1];
Q_a   = size(Ak, 1);

% FE- und RB-Systeme aufstellen
A_rb = Theta(1) * Ak_rb{1};
for k = 2:Q_a
    A_rb = A_rb + Theta(k) * Ak_rb{k};
end
A_fe = assemble_fe_A(Ak, Theta);

% Beide Lösungen bestimmen
U_rb           = A_rb \ f_rb;
U_rb_dirichlet = assemble_sol(Z * U_rb);
U_rb_sol       = Z * U_rb;
U_fe           = A_fe \ f;
U_fe_dirichlet = assemble_sol(U_fe);

% Output-Funktional für beide bestimmen und Differenz bestimmen
S_rb   = f_rb' * U_rb;
S_fe   = f' * U_fe;
diff_S = S_fe - S_rb;

% Fehlerverifikation
[err_s_mat, err_en_mat, alpha_LB_mat] = estimate_errors(U_rb, Theta, G);

disp(['Parameter           mu = ', mat2str(online_mu)]);
disp(['FE-Lösung:       s(mu) = ', num2str(S_fe, '%10.6f')]);
disp(['RB-Lösung:     s_N(mu) = ', num2str(S_rb, '%10.6f')]);
disp(['Echter Fehler:           ', num2str(diff_S, '%10.6e\n')]);
disp(['Fehlerschätzer:          ', num2str(err_s_mat, '%10.6e\n')]);

if tgl_plot == 1
    figure
    pdesurf(p, t, U_rb_dirichlet);
    title('RB-Approximation')
    xlabel('x');
    ylabel('y');
    colormap jet;

    figure
    pdesurf(p, t, U_fe_dirichlet);
    title('FE-Approximation')
    xlabel('x');
    ylabel('y');
    colormap jet;

    figure
    pdesurf(p, t, abs(U_rb_dirichlet - U_fe_dirichlet));
    title('Differenz FE-RB')
    xlabel('x');
    ylabel('y');
    colormap jet;
end
