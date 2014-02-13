%% estimate_error_prepare: Berechnet den parameterunabhängigen Teil des
%% Fehlerschätzers.
function G = prepare_error_estimator(X, Ak, F, Z)
    % Größen
    Q_a = size(Ak, 1);
    N_fe = size(Ak{1}, 1);
    N_rb = size(Z, 2);

    % parameterunabhänginger Anteil des Fehlerschätzers bestimmen
    H = zeros(N_fe, 1 + Q_a * N_rb);
    H(:, 1) = F;
    for k = 1:Q_a
        H(:, (1+k):Q_a:end) = - Ak{k} * Z;
    end
    G = H.' * (X.' \ H);
end
