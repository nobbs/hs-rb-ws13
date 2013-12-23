function [U] = rb_online(Mu, Z, Ak, F)
A_rb = zeros(size(Z, 2));
P = length(Mu);
F_rb = Z' * F;

for i = 1:P
	A_rb = A_rb + Mu(i) * Z' * Ak(:, :, i) * Z;
end

A_rb = A_rb + Z' * Ak(:, :, P + 1) * Z;
U_rb = A_rb \ F_rb
U = Z * U_rb;
end
