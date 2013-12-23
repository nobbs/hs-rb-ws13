% Liefert die exakte Lösung für Grid und eps
function [U] = exakteLoesung(T, ep)
	U = (- ( - exp(1 / ep) * T + exp(T / ep) + T - 1) / (exp(1 / ep) - 1))';
end
