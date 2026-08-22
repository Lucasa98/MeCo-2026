function x = coeficientes_descentrada(k,p,izq = false)
  % k: orden de la derivada
  % p: orden de precision

  % numero de puntos (que no sean el central)
  N = k+p-1;
  nodos = [1:N];

  % Matriz de coeficientes
  M = zeros(N, N);

  for i = 1:N
    M(i,:) = nodos.^i;
  endfor

  b = zeros(N, 1);
  b(k) = factorial(k);

  % Resolver sistema y agregar coeficiente de nodo 0
  c = M \ b;
  x = [-sum(c); c(:)];

  % descentrada izquierda: cambia de signo si k es impar
  if(izq)
    x = (-1)^k * flip(x);
  endif
endfunction
