X_inv = X^-1;

Q_a = size(Ak, 1);

AkF = cell(Q_a, N);

for m = 1:Q_a
    for n = 1:N
        AkF{m, n} = F' * X_inv * (Ak{m} * Z(:, n));
    end
end

AkAk = cell(Q_a, Q_a, N, N);
for m = 1:Q_a
    for q = 1:Q_a
        for n = 1:N
            for k = 1:N
                AkAk{m, q, n, k} = ((Ak{m} * Z(:, n))' * X_inv * (Ak{q} * Z(:, k)));
            end
        end
    end
end
