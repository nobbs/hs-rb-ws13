%% test_errors: function description
function [Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, F, X, Z, Xi_test, mu_bar, choosen_mus)
	global tgl_plot;
	global tgl_test;
	global tgl_pause;
	global tgl_logplot;
	global tgl_print;

	P = length(Ak);
	Q_a = P;
	N = size(Z, 2);
	num_ps = size(Ak{1}, 1);
	n_test = size(Xi_test, 2);

	% Immer wieder benötigt...
	Ak_rb = cell(P, 1);
	for j = 1:P
		Ak_rb{j} = Z' * Ak{j} * Z;
	end
	F_rb = Z' * F;

	H = zeros(num_ps, 1 + Q_a * N);
	H(:, 1) = F;
	for k = 1:Q_a
		H(:, (1+k):Q_a:end) = - Ak{k} * Z;
    end
    H = sparse(H);
	G = H.' * (X.' \ H);

	Test_errors_s_N = zeros(n_test, 1);
	Test_S_diff = zeros(n_test, 1);
	Test_effectivities = zeros(n_test, 1);

	if tgl_plot == 1
		Test_errors_u_en_N = zeros(n_test, 1);
		Test_u_en_diff = zeros(n_test, 1);
		Test_errors_u_X_N = zeros(n_test, 1);
		Test_u_X_diff = zeros(n_test, 1);
	end

	for j = 1:n_test
		if mod(j, 1000) == 0
			fprintf('-- Fehler bestimmt für %6d / %6d Test-Parameter\n', j, n_test);
		end

		mu = Xi_test(:, j);
		Theta = [mu; 1];

		A_rb = Theta(1) * Ak_rb{1};
		for k = 2:P
			A_rb = A_rb + Theta(k) * Ak_rb{k};
		end
		U_rb = A_rb \ F_rb;

		%%% Fehler berechnen
		Eps = zeros(1 + Q_a * N, 1);
		Eps(1) = 1;
		for k = 1:Q_a
			Eps((1+k):Q_a:end) = Theta(k) * U_rb;
        end
        Eps = sparse(Eps);
		residue_ = sqrt(abs(Eps.' * G * Eps));

		alpha_LB = min([mu ./ mu_bar; 1]);
		err_en = residue_ / sqrt(alpha_LB);
		err_s = (Eps.' * G * Eps) / alpha_LB;

		Test_errors_s_N(j) = err_s;

		% Effektivität
		A_fe = fe_assemble_A(Ak, mu);
		U_fe = A_fe \ F;

		S_fe = F.' * U_fe;
		S_rb = F_rb.' * U_rb;

		Test_effectivities(j) = err_s / (S_fe - S_rb);
		Test_S_diff(j) = (S_fe - S_rb);

		if tgl_plot == 1
			Test_errors_u_en_N(j) = err_en;
			Test_errors_u_X_N(j) = residue_ / alpha_LB;
			diff_u = (Z * U_rb - U_fe);
			Test_u_en_diff(j) = diff_u' * A_fe * diff_u;
			Test_u_X_diff(j) = diff_u' * X * diff_u;
		end
	end

	% Test_errors_s_N
	% Test_effectivities
	% Test_S_diff

	if tgl_plot == 1
		figure()
		% subplot(1, 2, 1)
		if tgl_logplot == 1
			semilogy(Xi_test(1:10:end), [abs(Test_errors_s_N(1:10:end)), abs(Test_S_diff(1:10:end))]);
		else
			plot(Xi_test(1:10:end), [abs(Test_errors_s_N(1:10:end)), abs(Test_S_diff(1:10:end))]);
		end
		if tgl_print == 1
			children = get(gca, 'Children');
			set(children(1), 'Color', 'b', 'LineStyle', '--');
			set(children(2), 'Color', 'r');
		else
			title('Fehler des Funktionals s(\mu)');
			ylabel('Abweichung');
			xlabel('Parameter \mu');
			% ylim([1e-16 1e0]);
			legend('Approx. Fehler', 'Tatsächlicher Fehler', 'Location', 'SouthEast');
		end

		if tgl_print == 0
			figure()
			% subplot(1, 2, 2)
			if tgl_logplot == 1
				semilogy(Xi_test, [Test_errors_u_en_N, abs(Test_u_en_diff), abs(Test_errors_u_X_N), abs(Test_u_X_diff)]);
			else
				plot(Xi_test, [Test_errors_u_en_N, abs(Test_u_en_diff), abs(Test_errors_u_X_N), abs(Test_u_X_diff)]);
			end

			if tgl_print == 1
				children = get(gca, 'Children');
				set(children(1), 'Color', 'c', 'LineStyle', '--');
				set(children(2), 'Color', 'g', 'LineStyle', '-.');
				set(children(3), 'Color', 'b', 'LineStyle', ':');
				set(children(4), 'Color', 'r');
			else
				children = get(gca, 'Children');
				set(children(1), 'Color', 'r', 'LineStyle', '--');
				set(children(2), 'Color', 'r');
				set(children(3), 'Color', 'b', 'LineStyle', '--');
				set(children(4), 'Color', 'b');
			end

			title('Fehler der Lösung u(\mu)');
			ylabel('Abweichung');
			xlabel('Parameter \mu');
			% ylim([1e-16 1e0]);
			legend('Approx. Fehler in Energie-Norm', 'Tatsächlicher Fehler in Energie-Norm', 'Approx. Fehler in X-Norm', 'Tatsächlicher Fehler in X-Norm', 'Location', 'SouthEast');
		end

		if tgl_pause == 1
			pause('on');
			pause;
		else
			pause('off');
			pause(2);
		end
	end

	% figure();
	% semilogy(Xi_test, [Test_errors_s_N, abs(Test_S_diff), Test_errors_u_en_N, abs(Test_u_en_diff)]);
	% title('Approximiertier und tatsächlicher Fehler des Funktionals / der Lösung');
	% ylabel('Abweichung');
	% xlabel('Parameter');
	% ylim([1e-16 1e0]);
	% legend('Approx. Fehler', 'Tatsächlicher Fehler', 'Approx. Fehler', 'Tatsächlicher Fehler');
	% waitforbuttonpress;

	% 2 Parameter
	% figure()
	% n_reshape = floor(sqrt(n_test));
	% X_mesh = reshape(Xi_test(1, :), n_reshape, n_reshape);
	% Y_mesh = reshape(Xi_test(2, :), n_reshape, n_reshape);
	% surf(X_mesh, Y_mesh, reshape(Test_errors_s_N, n_reshape, n_reshape));
	% hold on;
	% surf(X_mesh, Y_mesh, reshape(abs(Test_S_diff), n_reshape, n_reshape));
	% hold off;
	% pause(1);

	Delta_s_N_max = max(Test_errors_s_N);
	eta_s_N_ave = nanmean(Test_effectivities(isfinite(Test_effectivities)));
	eta_s_N_max = max(Test_effectivities(isfinite(Test_effectivities)));
	rho_S_err_N = Delta_s_N_max / max(Test_S_diff);
end
