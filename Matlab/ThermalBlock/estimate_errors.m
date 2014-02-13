%% estimate_errors: Berechnet die drei Fehlerschätzer und alpha_LB für die
%% gegebene RB-Lösung und den Parameter.
function [err_s, err_mu, err_X, alpha_LB] = estimate_errors(u_rb, Theta_a, G)
    % Größen
    Q_a  = size(Theta_a, 1);
    N_rb = size(u_rb, 1);

    % parameterabhängigen Teil des Fehlerschätzers aufbauen
    Eps    = zeros(1 + Q_a * N_rb, 1);
    Eps(1) = 1;
    for k = 1:Q_a
        Eps((1+k):Q_a:end) = Theta_a(k) * u_rb;
    end

    e_hat_X  = sqrt(abs(Eps.' * G * Eps));
    alpha_LB = min(Theta_a);

    % Fehlerschätzer ausrechnen
    err_mu = e_hat_X / sqrt(alpha_LB);
    err_X  = e_hat_X / alpha_LB;
    err_s  = e_hat_X^2 / alpha_LB;
end
