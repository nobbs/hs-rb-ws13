%% test_errors: Bestimmt den FE-RB-Fehler mittels Sch�tzer und L�sen mittels
%% beider Verfahren. Plotten / Ausgeben der Fehlersch�tzer, Effiktiv�ten.
function [Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, f, X, Z, Xi_test, mu_bar)
	% Globale Flags / Optionen
	global tgl_plot;
	global tgl_test;
	global tgl_pause;
	global tgl_logplot;

	% Gr��en
	Q_a    = length(Ak);
	N_fe   = size(Z, 1);
	N_rb   = size(Z, 2);
	n_test = size(Xi_test, 2);

	% Parameterunabh�ngiger Anteil des RB-Gleichungssystems
	Ak_rb = cell(Q_a, 1);
	for j = 1:Q_a
		Ak_rb{j} = Z' * Ak{j} * Z;
	end
	f_rb = Z' * f;

	% Parameterunabh�ngiger Fehlersch�tzer-Anteil
	G = prepare_error_estimator(X, Ak, f, Z);

	% Vektoren um Fehler f�r jeden Testparameter abzuspeichern.
	Stat_err_s  = zeros(n_test, 1);
	Stat_diff_s = zeros(n_test, 1);
	Stat_eff_s  = zeros(n_test, 1);

	% Wenn geplottet werden soll, dann werden noch mehr Fehler bestimmt.
	if tgl_plot == 1
		Stat_err_u_mu  = zeros(n_test, 1);
		Stat_diff_u_mu = zeros(n_test, 1);
		Stat_err_u_X   = zeros(n_test, 1);
		Stat_diff_u_X  = zeros(n_test, 1);
	end

	for j = 1:n_test
		% Progress-Ausgabe
        if mod(j, n_test / 10) == 0
			fprintf('-- Fehler bestimmt f�r %6d / %6d Test-Parameter\n', j, n_test);
		end

		% Testparameter mu w�hlen und Theta(mu) bestimmen
		mu    = Xi_test(:, j);
		Theta = [mu; 1];

		% RB-System zusammensetzen und RB-L�sung bestimmen
		A_rb = Theta(1) * Ak_rb{1};
		for k = 2:Q_a
			A_rb = A_rb + Theta(k) * Ak_rb{k};
		end
		U_rb = A_rb \ f_rb;

		% Fehlersch�tzer berechnen
		[err_s, err_mu, err_X, alpha_LB] = estimate_errors(U_rb, Theta, G);

		% Finite-Elemente-L�sung bestimmen
		A_fe = assemble_fe_A(Ak, [mu, 1]);
		U_fe = A_fe \ f;

		% Output-Funktionale f�r RB und FE auswerten
		s_rb = f_rb.' * U_rb;
		s_fe = f.' * U_fe;

		% Fehler und Fehlersch�tzer abspeichern
		Stat_err_s(j)  = err_s;
		Stat_eff_s(j)  = err_s / (s_fe - s_rb);
		Stat_diff_s(j) = (s_fe - s_rb);

		if tgl_plot == 1
			diff_u            = (Z * U_rb - U_fe);
			Stat_err_u_mu(j)  = err_mu;
			Stat_err_u_X(j)   = err_X;
			Stat_diff_u_mu(j) = diff_u' * A_fe * diff_u;
			Stat_diff_u_X(j)  = diff_u' * X * diff_u;
		end
	end

	% Plotten, falls erw�nscht
	if tgl_plot == 1
		figure
		% subplot(1, 2, 1)

		% Logarithmisch oder linear?
		if tgl_logplot == 1
			semilogy(Xi_test(1:10:end), [abs(Stat_err_s(1:10:end)), abs(Stat_diff_s(1:10:end))]);
		else
			plot(Xi_test(1:10:end), [abs(Stat_err_s(1:10:end)), abs(Stat_diff_s(1:10:end))]);
		end

		title('Fehler des Funktionals s(\mu)');
		ylabel('Abweichung');
		xlabel('Parameter \mu');
		legend('Approx. Fehler', 'Tats�chlicher Fehler', 'Location', 'SouthEast');

		figure
		% subplot(1, 2, 2)

		% Logarithmisch oder linear?
		if tgl_logplot == 1
			semilogy(Xi_test, [Stat_err_u_mu, abs(Stat_diff_u_mu), abs(Stat_err_u_X), abs(Stat_diff_u_X)]);
		else
			plot(Xi_test, [Stat_err_u_mu, abs(Stat_diff_u_mu), abs(Stat_err_u_X), abs(Stat_diff_u_X)]);
		end

		children = get(gca, 'Children');
		set(children(1), 'Color', 'r', 'LineStyle', '--');
		set(children(2), 'Color', 'r');
		set(children(3), 'Color', 'b', 'LineStyle', '--');
		set(children(4), 'Color', 'b');

		title('Fehler der L�sung u(\mu)');
		ylabel('Abweichung');
		xlabel('Parameter \mu');
		legend('Approx. Fehler in \mu-Norm', 'Tats�chlicher Fehler in \mu-Norm', 'Approx. Fehler in X-Norm', 'Tats�chlicher Fehler in X-Norm', 'Location', 'SouthEast');

		if tgl_pause == 1
			pause('on');
			pause;
		else
			pause('off');
			pause(2);
		end
	end

	% Maximalen Fehler und Effektivit�ten bestimmen (und dabei NaN / inf
	% ignorieren).
	Delta_s_N_max = max(Stat_err_s);
	eta_s_N_ave   = nanmean(Stat_eff_s(isfinite(Stat_eff_s)));
	eta_s_N_max   = max(Stat_eff_s(isfinite(Stat_eff_s)));
	rho_S_err_N   = Delta_s_N_max / max(Stat_diff_s);
end
