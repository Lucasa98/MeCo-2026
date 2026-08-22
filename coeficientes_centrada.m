function x = coeficientes_centrada(k, p)
  % k: orden de la derivada >= 1
  % p: orden de la aproximacion >= 1

  % numero de puntos (que no sean el central)
  N = max(2, ceil((p + k - 2) / 2) * 2);

  % Matriz de coeficientes
  M = ones(k+1, N);
  b = zeros(k+1, 1);
  b(k) = factorial(k);

  for i = 1:k+1
    for j = 1:(N/2)
      M(i,2*j-1) = factorial(j)^i;
      M(i,2*j) = (-factorial(j))^i;
    endfor
  endfor

  x = M \ b;
endfunction
