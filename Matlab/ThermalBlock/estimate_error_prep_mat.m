X_inv = X^-1;

Q_a = size(Ak, 1);
num_ps = size(Ak{1}, 1);

H = zeros(num_ps, 1 + Q_a * N);
H(:, 1) = F;
for k = 1:Q_a
    H(:, (1+k):Q_a:end) = - Ak{k} * Z;
end
H = sparse(H);
G = H.' * (X.' \ H);