function [err_s, err_en, alpha_LB] = estimate_error_mat(u_rb, Theta_a, G, Q_a, N)
%     Q_a = size(AkF, 1);
%     N = size(AkF, 2);

    Eps = zeros(1 + Q_a * N, 1);
    Eps(1) = 1;
    for k = 1:Q_a
        Eps((1+k):Q_a:end) = Theta_a(k) * u_rb;
    end
    Eps = sparse(Eps);
    residue_ = abs(Eps.' * G * Eps);

%     err_ = err_ff - 2 * err_fa + err_aa;
    err_ = residue_;
    alpha_LB = min(Theta_a);

    err_en = sqrt(abs(err_ / alpha_LB));
    err_s = err_en^2;
end
