function x = coeficientes_centrada(k, p)
  % coeficientes para armar el stencil de -m*h hasta +m*h
  % k: orden de la derivada >= 1
  % p: orden de la aproximacion >= 1

  % numero de puntos (que no sean el central)
  N = 2*floor((k+1)/2) + 2*ceil(p/2) - 2;
  m = N/2;
  nodos = [-m:-1, 1:m];

  % Matriz de coeficientes
  M = zeros(N, N);

  for i = 1:N
    M(i,:) = nodos.^i;
  endfor

  b = zeros(N, 1);
  b(k) = factorial(k);

  c = M \ b;
  x = [c(1:m); -sum(c); c(m+1:end)];
endfunction
