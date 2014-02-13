% Statistik-Speicher
tab_N             = [];
tab_Delta_s_N_max = [];
tab_eta_s_N_ave   = [];
tab_eta_s_N_max   = [];
tab_rho_S_err_N   = [];

% Größen
N    = 1;
N_fe = size(Ak{1}, 1);
Q_a  = size(Ak, 1);

% Speicher für gewählte Parameter (Index bzw. Parameter mu bzw. u(mu)).
snapshots_idx = [];
snapshots     = sparse(N_fe, 1);
Z             = sparse(N_fe, 1);

% Ersten Snapshot berechnen (für zufällig gewählten Parameter bzw. mu = mu_bar)
mu_ss = mu_bar;

% idx_ss = randi(n_train, 1);
% mu_ss = Xi_train(:, idx_ss);
% snapshots_idx(1) = idx_ss;

% Gram-Schmidt-Verfahren für bessere Kondition
snapshots(:, 1) = assemble_fe_A(Ak, [mu_ss, 1]) \ f;
Z(:, 1)         = snapshots(:, 1) / sqrt(snapshots(:, 1)' * X * snapshots(:, 1));

disp(['<< N = ', mat2str(N), ', mu = ', mat2str(mu_ss)])

% Fehlerschätzer mit tatsächlichem FE-RB-Fehler vergleichen
if tgl_test ~= 0
    disp('<< Testergebnisse')
    [Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, f, X, Z, Xi_test, mu_bar);

    tab_N                      = [N];
    tab_Delta_s_N_max(end + 1) = Delta_s_N_max;
    tab_eta_s_N_ave(end + 1)   = eta_s_N_ave;
    tab_eta_s_N_max(end + 1)   = eta_s_N_max;
    tab_rho_S_err_N(end + 1)   = rho_S_err_N;

    print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);
end

% Greedy-Verfahren starten
disp(['<< Prüfe N + 1 = ', num2str(N + 1)])
for i = 2:N_max
    disp('<< Bestimme nächsten Snapshot...')

    idx_ss = -1;        % Index des neuen Trainings-Parameter
    mu_ss  = -1;        % neuer Trainings-Parameter
    err_mu = 0;         % geschätzter Fehler in mu-Norm
    err_X  = 0;         % geschätzter Fehler in X-Norm
    err_s  = 0;         % geschätzter Funktional-Fehler
    u_N    = 0;         % RB-Lösung
    s_N    = 0;         % RB-Funktionalwert

    % Parameterunabhängiger Fehlerschätzer-Anteil
    G = prepare_error_estimator(X, Ak, f, Z);

    % Parameterunabhängiger Anteil des RB-Gleichungssystems
    Ak_rb = cell(Q_a, 1);
    for j = 1:Q_a
        Ak_rb{j} = Z' * Ak{j} * Z;
    end
    f_rb = Z' * f;

    % Suche maximalen Fehler über alle Trainings-Parameter
    for j = 1 : size(Xi_train, 2)
        % Ignoriere bereits verwendete Parameter
        if any(snapshots_idx == j)
            continue;
        end

        % Progress-Ausgabe
        if mod(j, n_train / 10) == 0
            fprintf('-- Fehler bestimmt für %6d / %6d Trainings-Parameter\n', j, n_train);
        end

        % Trainings-Parameter holen
        mu_check = Xi_train(:, j);
        Theta    = [mu_check; 1];

        % RB-System zusammenstellen und Lösung bestimmen
        A_rb = Ak_rb{Q_a};
        for k = 1:Q_a-1
            A_rb = A_rb + Theta(k) * Ak_rb{k};
        end
        u_rb = A_rb \ f_rb;

        % Für Fehlerschätzer benötigt
        [err_s_check, err_mu_check, err_X_check, alpha_LB] = estimate_errors(u_rb, Theta, G);

        % Prüfen, ob berechneter Fehler größer als bisheriges Maximum
        if err_mu_check > err_mu
            % RB-Lösung auf passendes Format bringen und Output-Funktional auswerten
            u_N = Z * u_rb;
            s_N = f_rb' * u_rb;

            % Neues mu mit max. Fehlerschätzer abspeichern
            idx_ss = j;
            mu_ss  = mu_check;

            % Fehler abspeichern
            err_mu = err_mu_check;
            err_s  = err_s_check;
            err_X  = err_X_check;
        end
    end

    % Prüfen, ob die gewünschte Toleranz bereits erreicht wurde.
    if err_s < tolerance
        fprintf('Vorgegebene Fehlertoleranz %10.6e für s bereits erreicht mit N = %d\n', tolerance, N);
        fprintf('Delta^en_N_max =  %8.6e\n', full(err_mu))
        fprintf('Delta_N_max    =  %8.6e\n', full(err_X))
        fprintf('Delta^s_N_max  =  %8.6e\n\n', full(err_s))
        print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);
        break;
    end

    fprintf(['\nSnapshot gewählt mit mu = ', mat2str(mu_ss)])

    % FE-Lösung zu neu gewähltem mu bestimmen und Gram-Schmidt darauf anwenden
    A_ss = assemble_fe_A(Ak, [mu_ss, 1]);
    snapshots(:, i)        = sparse(A_ss \ f);
    snapshots_idx(end + 1) = idx_ss;

    z = snapshots(:, i);
    for j = 1:(i-1)
        z = z - (snapshots(:, i)' * X * Z(:, j)) * Z(:, j);
    end
    Z(:, i) = z / (sqrt(z' * X * z));

    % Noch ein paar Ausgaben
    diff_u_mu = sqrt(full((snapshots(:, i) - u_N).' * A_ss * (snapshots(:, i) - u_N)));
    diff_u_X  = sqrt(full((snapshots(:, i) - u_N).' * X * (snapshots(:, i) - u_N)));
    s_fe      = full(f.' * snapshots(:, i));

    fprintf('\n\n||u - u_N||_mu  =  %8.6e\n', full(diff_u_mu))
    fprintf('Delta^en_N_max  =  %8.6e\n\n', full(err_mu))
    fprintf('||u - u_N||_X   =  %8.6e\n', full(diff_u_X))
    fprintf('Delta_N_max     =  %8.6e\n\n', full(err_X))
    fprintf('  s - s_N       =  %8.6e\n', full(s_fe - s_N))
    fprintf('Delta^s_N_max   =  %8.6e\n\n', full(err_s))

    N = i;
    disp(['<< Prüfe N + 1 = ', num2str(N + 1)])

    % Noch mehr Ausgaben
    if tgl_test ~= 0
        disp('<< Testergebnisse')
        [Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, f, X, Z, Xi_test, mu_bar);

        tab_N(end + 1)             = [N];
        tab_Delta_s_N_max(end + 1) = Delta_s_N_max;
        tab_eta_s_N_ave(end + 1)   = eta_s_N_ave;
        tab_eta_s_N_max(end + 1)   = eta_s_N_max;
        tab_rho_S_err_N(end + 1)   = rho_S_err_N;

        print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);
    end
end

% Tabelle als Cell-Array um damit eine LaTeX-Tabelle generieren zu können
Table = {tab_N', tab_Delta_s_N_max', tab_eta_s_N_ave', tab_eta_s_N_max', tab_rho_S_err_N'};
