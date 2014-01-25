% Um die Tabelle zu erzeugen
tab_N = [];
tab_Delta_s_N_max = [];
tab_eta_s_N_ave = [];
tab_eta_s_N_max = [];
tab_rho_S_err_N = [];

% Immer wieder benötigt...
Choosen_idx = [];
Snapshots = sparse(num_ps_no_dbc, 1);
Z = sparse(num_ps_no_dbc, 1);
num_Ak = size(Ak, 1);

% Ersten Snapshot berechnen (für zufällig gewählten Parameter)
% tic;
idx_ss = randi(n_train, 1);
% mu_ss = Xi_train(:, idx_ss);
% Choosen_idx(1) = idx_ss;
mu_ss = mu_bar;
% und mittels Gram-Schmidt Orthonormalisieren (-> besser konditioniert...)
Snapshots(:, 1) = fe_assemble_A(Ak, mu_ss) \ F;
Z(:, 1) = Snapshots(:, 1) / sqrt(Snapshots(:, 1)' * X * Snapshots(:, 1));

N = 1;
disp(['<< N = ', mat2str(N), ', mu = ', mat2str(mu_ss)])
% toc
disp('<< Testergebnisse')
if tgl_test ~= 0
	% tic
	[Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, F, X, Z, Xi_test, mu_bar, Xi_train(Choosen_idx));
	tab_N = [1];
	tab_Delta_s_N_max(end + 1) = Delta_s_N_max;
	tab_eta_s_N_ave(end + 1) = eta_s_N_ave;
	tab_eta_s_N_max(end + 1) = eta_s_N_max;
	tab_rho_S_err_N(end + 1) = rho_S_err_N;
	print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);
	% toc
end
disp([' ']);

disp(['<< Prüfe N + 1 = ', num2str(N + 1)])
for i = 2:N_max
	% tic
	disp('<< Bestimme nächsten Snapshot...')

	idx_ss = -1;		% Index des neuen Trainings-Parameter
	mu_ss = -1;			% neuer Trainings-Parameter
	err_en = 0;			% geschätzter Fehler in Energie-Norm
	err_X = 0;			% geschätzter Fehler in X-Norm
	err_s = 0;			% geschätzter Funktional-Fehler
	U_rb_N = 0;			% RB-Lösung
	S_N = 0;			% RB-Funktionalwert

	% Benötigt für Fehlerschätzer
	H = zeros(num_ps_no_dbc, 1 + (bx * by) * N);
	H(:, 1) = F;
	for k = 1:num_Ak
		H(:, (1+k):num_Ak:end) = - Ak{k} * Z;
    end
    H = sparse(H);
	G = H.' * (X.' \ H);

	% Benötigt für RB-FEM
	Ak_rb = cell(num_Ak, 1);
	for j = 1:num_Ak
		Ak_rb{j} = Z' * Ak{j} * Z;
	end
	F_rb_check = Z' * F;

	% Speichere geschätzte Fehler für Plot
	Errors_en = zeros(size(Xi_train, 2), 1);

	% suche maximalen Fehler über alle Trainings-Parameter
	for j = 1 : size(Xi_train, 2)
		% Ignoriere bereits verwendete Parameter
		if any(Choosen_idx == j)
			continue;
		end

		if mod(j, n_train / 10) == 0
			fprintf('-- Fehler bestimmt für %6d / %6d Trainings-Parameter\n', j, n_train);
		end

		% Trainings-Parameter
		mu_check = Xi_train(:, j);
		Theta = [mu_check; 1];

		% RB-Lösung bestimmen
		A_rb_check = Ak_rb{num_Ak};
		for k = 1:num_Ak-1
			A_rb_check = A_rb_check + Theta(k) * Ak_rb{k};
		end
		U_rb_check = A_rb_check \ F_rb_check;

		% Für Fehlerschätzer benötigt
		Eps = zeros(1 + (bx * by) * N, 1);
		Eps(1) = 1;
		for k = 1:num_Ak
			Eps((1+k):num_Ak:end) = Theta(k) * U_rb_check;
        end
        Eps = sparse(Eps);

		% Fehler schätzen
		err_residue = sqrt(abs(Eps.' * G * Eps));
		alpha_LB = min([mu_check ./ mu_bar; 1]);

		err_en_check = err_residue / sqrt(alpha_LB);
		err_s_check = err_en_check^2;

		Errors_en(j) = err_en_check;

		if err_en_check > err_en
			U_rb_N = Z * U_rb_check;
			S_N = F_rb_check' * U_rb_check;

			idx_ss = j;
			mu_ss = mu_check;

			err_en = err_en_check;
			err_s = err_s_check;
			err_X = err_en_check / sqrt(alpha_LB);
		end
	end

	if err_s < tolerance
		fprintf('Vorgegebene Fehlertoleranz %10.6e für s bereits erreicht mit N = %d\n', tolerance, N);
		disp(' ')
		% fprintf('||u - u_N||_mu  =  %8.6e\n', full(diff_u_en_norm))
		fprintf('Delta^en_N_max  =  %8.6e\n\n', full(err_en))
		% fprintf('||u - u_N||_X   =  %8.6e\n', full(diff_u_X_norm))
		fprintf('Delta_N_max     =  %8.6e\n\n', full(err_X))
		% fprintf('  s - s_N       =  %8.6e\n', full(S - S_N))
		fprintf('Delta^s_N_max   =  %8.6e\n', full(err_s))
		disp(' ')

		print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);
		break;
	end

	disp(['Snapshot gewählt mit mu = ', mat2str(mu_ss)])
	% A_ss = mu_ss * Ak(:, :, 1) + Ak(:, :, 2);
	A_ss = fe_assemble_A(Ak, mu_ss);
	Snapshots(:, i) = sparse(A_ss \ F);
	Choosen_idx(end + 1) = idx_ss;

	% Orthonormalisieren
	z = Snapshots(:, i);
	for j = 1:(i-1)
		z = z - (Snapshots(:, i)' * X * Z(:, j)) * Z(:, j);
	end
	Z(:, i) = z / (sqrt(z' * X * z));



	% aktuellen Stand ausgeben
	% disp('Fehlerschaetzer err_en, echter Fehler ||u - u_N||_mu und Effektivität eta^en_N')

	diff_u_en_norm = sqrt(full((Snapshots(:, i) - U_rb_N).' * A_ss * (Snapshots(:, i) - U_rb_N)));
	[err_en, diff_u_en_norm, err_en / diff_u_en_norm];

	% disp('Fehlerschaetzer err_X, echter Fehler ||u - u_N||_X und Effektivität eta_N')
	diff_u_X_norm = sqrt(full((Snapshots(:, i) - U_rb_N).' * X * (Snapshots(:, i) - U_rb_N)));
	[err_X, diff_u_X_norm, err_X / diff_u_X_norm];

	% disp('Fehlerschaetzer err_s, echter Fehler |s - s_N| und Effektivität eta^s_N')
	S = full(F.' * Snapshots(:, i));
	[err_s, S - S_N, err_s / (S - S_N)];

	disp(' ')
	fprintf('||u - u_N||_mu  =  %8.6e\n', full(diff_u_en_norm))
	fprintf('Delta^en_N_max  =  %8.6e\n\n', full(err_en))
	fprintf('||u - u_N||_X   =  %8.6e\n', full(diff_u_X_norm))
	fprintf('Delta_N_max     =  %8.6e\n\n', full(err_X))
	fprintf('  s - s_N       =  %8.6e\n', full(S - S_N))
	fprintf('Delta^s_N_max   =  %8.6e\n', full(err_s))
	disp(' ')

	% aktuellen Stand plotten
	% 2 Parameter
	% n_reshape = floor(sqrt(n_train));
	% surf(reshape(Xi_train(1, :), n_reshape, n_reshape), reshape(Xi_train(2, :), n_reshape, n_reshape), reshape(Errors_en, n_reshape, n_reshape));
	% hold on;
	% for y = Xi_train(:, Choosen_idx)
	% 	plot(y, '*r');
	% end
	% hold off;

	% 1 Paramater
	% plot(Xi_train, Errors_en);
	% hold on;
	% for y = Xi_train(:, Choosen_idx)
	% 	plot(y, 0, '*r');
	% end
	% hold off;
	% pause(1)

	% toc

	N = i;
	disp(['<< Prüfe N + 1 = ', num2str(N + 1)])

	if tgl_test ~= 0
		% tic
		[Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, F, X, Z, Xi_test, mu_bar, Xi_train(Choosen_idx));
		tab_N(end + 1) = N;
		tab_Delta_s_N_max(end + 1) = Delta_s_N_max;
		tab_eta_s_N_ave(end + 1) = eta_s_N_ave;
		tab_eta_s_N_max(end + 1) = eta_s_N_max;
		tab_rho_S_err_N(end + 1) = rho_S_err_N;
		print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);
		% toc
	end
end

Table = {tab_N', tab_Delta_s_N_max', tab_eta_s_N_ave', tab_eta_s_N_max', tab_rho_S_err_N'};
