% Immer wieder benötigt...
Choosen_idx = [];
Snapshots = sparse(num_ps_no_dbc, 1);
Z = sparse(num_ps_no_dbc, 1);
num_Ak = length(Ak);

% Ersten Snapshot berechnen (für zufällig gewählten Parameter)
idx_ss = randi(n_train, 1);
mu_ss = Xi_train(idx_ss);
Choosen_idx(1) = idx_ss;
% und mittels Gram-Schmidt Orthonormalisieren (-> besser konditioniert...)
Snapshots(:, 1) = fe_assemble_A(Ak, mu_ss) \ F;
Z(:, 1) = Snapshots(:, 1) / sqrt(Snapshots(:, 1)' * X * Snapshots(:, 1));

N = 1;
disp(['<< N = ', num2str(N), ', zufällig gewähltes mu = ', num2str(mu_ss)])
disp('<< Testergebnisse')
[Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, F, X, Z, Xi_test, mu_bar)

for i = 2:N_max
	disp('<< Bestimme nächsten Snapshot...')

	idx_ss = -1;		% Index des neuen Trainings-Parameter
	mu_ss = -1;			% neuer Trainings-Parameter
	err_en = 0;			% geschätzter Fehler
	err_s = 0;			% geschätzter Funktional-Fehler
	U_rb_N = 0;			% RB-Lösung
	S_N = 0;			% RB-Funktionalwert

	% Benötigt für Fehlerschätzer
	H = sparse(num_ps_no_dbc, 1 + 2 * N);
	H(:, 1) = F;
	for k = 1:num_Ak
		H(:, (1+k):num_Ak:end) = - Ak{k} * Z;
	end
	G = H' * (X \ H);

	% Benötigt für RB-FEM
	F_rb_check = Z' * F;

	% Speichere geschätzte Fehler für Plot
	Errors_en = zeros(size(Xi_train, 2), 1);

	% suche maximalen Fehler über alle Trainings-Parameter
	for j = 1 : size(Xi_train, 2)
		% Ignoriere bereits verwendete Parameter
		if any(Choosen_idx == j)
			continue;
		end

		if mod(j, 1000) == 0
			disp(['Fehler bestimmt für ', num2str(j), ' von ', num2str(n_train), ' Trainings-Parameter'])
		end

		% Trainings-Parameter
		mu_check = Xi_train(:, j);
		Theta = [mu_check; 1];

		% RB-Lösung bestimmen
		A_rb_check = Z' * (fe_assemble_A(Ak, mu_check)) * Z;
		U_rb_check = A_rb_check \ F_rb_check;

		% Für Fehlerschätzer benötigt
		Eps = sparse(1 + 2 * N, 1);
		Eps(1) = 1;
		for k = 1:num_Ak
			Eps((1+k):num_Ak:end) = Theta(k) * U_rb_check;
		end

		% Fehler schätzen
		err_residue = full(real(sqrt(Eps.' * G * Eps)));
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
		end
	end

	disp(['Snapshot gewählt mit mu = ', num2str(mu_ss)])
	% A_ss = mu_ss * Ak(:, :, 1) + Ak(:, :, 2);
	A_ss = fe_assemble_A(Ak, mu_ss);
	Snapshots(:, i) = A_ss \ F;
	Choosen_idx(end + 1) = idx_ss;

	% Orthonormalisieren
	z = Snapshots(:, i);
	for j = 1:(i-1)
		z = z - (Snapshots(:, i)' * X * Z(:, j)) * Z(:, j);
	end
	Z(:, i) = z / (sqrt(z' * X * z));

	% aktuellen Stand ausgeben
	disp('Fehlerschaetzer err_en')
	err_en / sqrt(alpha_LB)
	disp('Echter Fehler ||u - u_N||_X')
	sqrt(full((Snapshots(:, i) - U_rb_N).' * A_ss * (Snapshots(:, i) - U_rb_N)))
	disp('Fehlerschaetzer err_s')
	err_s
	disp('Echter Fehler |s - s_N|')
	S = full(F.' * Snapshots(:, i));
	S - S_N
	disp('Effectivness eta^s_N')
	err_s / (S - S_N)

	% aktuellen Stand plotten
	plot(Xi_train, Errors_en);
	hold on;
	for y = Xi_train(Choosen_idx)
		plot(y, 0, '*r');
	end
	hold off;
	% xlim([mu_min, mu_max]);
	% ylim([1e-6, 100]);
	pause(1)
	% waitforbuttonpress;

	N = i;
	[Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, F, X, Z, Xi_test, mu_bar)

	if Delta_s_N_max < tolerance
		disp(['Vorgegebene Fehlertoleranz erreicht mit N = ', num2str(N), '.']);
		break;
	end
end
