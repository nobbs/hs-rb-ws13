function [err_s, err_en, alpha_LB] = estimate_error(u_rb, Theta_a, F, AkF, AkAk, X_inv)
    Q_a = size(AkF, 1);
    N = size(AkF, 2);

    % Fehlerschätzer
    err_ff = F' * X_inv * F;
    err_fa = 0;
    err_aa = 0;

    for m = 1:Q_a
        for n = 1:N
            err_fa = err_fa + Theta_a(m) * u_rb(n) * AkF{m, n};
        end
    end

    for m = 1:Q_a
        for q = 1:Q_a
            for n = 1:N
                for k = 1:N
                    err_aa = err_aa + Theta_a(m) * Theta_a(q) * u_rb(n) * u_rb(k) * AkAk{m, q, n, k};
                end
            end
        end
    end

    err_ = err_ff - 2 * err_fa + err_aa;
    alpha_LB = min(Theta_a);

    err_en = sqrt(abs(err_ / alpha_LB));
    err_s = err_en^2;
end
