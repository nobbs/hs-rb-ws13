% Parameter wählen
Theta = [online_mu; 1];
Q_a = size(Ak, 1);

% FE- und RB-Gls
A_rb = Ak_rb{num_Ak};
for i = 1:P
    A_rb = A_rb + online_mu(i) * Ak_rb{i};
end
A = fe_assemble_A(Ak, online_mu);

U_rb = A_rb \ F_rb;
U_rb_d = assemble_sol(Z * U_rb);
U_rb_sol = Z * U_rb;
U_fe = A \ F;
U_fe_d = assemble_sol(U_fe);

S_rb = F_rb' * U_rb;
S_fe = F' * U_fe;
err_true_s = S_fe - S_rb;

% Fehlerverifikation
if tgl_certify == 1
%     [err_s, err_en, alpha_LB] = estimate_error(U_rb, Theta, F, AkF, AkAk, X_inv);
    [err_s_mat, err_en_mat, alpha_LB_mat] = estimate_error_mat(U_rb, Theta, G, Q_a, N);
end

disp(['Parameter           mu = ']);
online_mu
disp(['FE-Lösung:       s(mu) = ', num2str(S_fe, '%10.6f')]);
disp(['RB-Lösung:     s_N(mu) = ', num2str(S_rb, '%10.6f')]);
disp(['Echter Fehler:           ', num2str(err_true_s, '%10.6e\n')]);
if tgl_certify == 1
%     disp(['Fehlerschätzer:          ', num2str(err_s, '%10.6e\n')]);
    disp(['Fehlerschätzer:          ', num2str(err_s_mat, '%10.6e\n')]);
end

if tgl_plot == 1
    figure
    pdesurf(p, t, U_rb_d);
    title('RB-Approximation')
    xlabel('x');
    ylabel('y');
    colormap jet;

    figure
    pdesurf(p, t, U_fe_d);
    title('FE-Approximation')
    xlabel('x');
    ylabel('y');
    colormap jet;

    figure
    pdesurf(p, t, abs(U_rb_d - U_fe_d));
    title('Differenz FE-RB')
    xlabel('x');
    ylabel('y');
    colormap jet;
end
