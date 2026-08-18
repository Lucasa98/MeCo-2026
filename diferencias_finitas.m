function [x,y] = diferencias_finitas(f, inter, y0, N)
  % diferencias finitas
  % f = [p(x); q(x); r(x)]
  % inter = [a b]
  % y0 = [y(a) y(b)]

  % h = (b-a)/(N+1)
  h = (inter(2) - inter(1))/(N+1);

  % construir sistema
  i = [1:N];
  xi = inter(1) + i*h;
  pqr = f(xi);
  a_i = 2 + h^2 * pqr(2,:);     % diagonal principal
  b_i = -1 + (h/2) * pqr(1,:);  % diagonal derecha
  c_i = -1 - (h/2) * pqr(1,:);  % diagonal izquierda
  d_i = -h^2 * pqr(3,:);

  % primer y ultima ecuacion
  d_i(1) = d_i(1) + (1 + (h/2) * pqr(1,1))*y0(1);
  d_i(N) = d_i(N) + (1 - (h/2) * pqr(1,N))*y0(2);

  % matriz
  M = diag(a_i);                % diagonal principal
  M = M + diag(c_i(2:N), -1);   % diagonal izquierda
  M = M + diag(b_i(1:N-1), 1);  % diagonal derecha

  % resolver sistema
  w = M \ d_i(:);

  % solucion con extremos
  x = [inter(1), xi, inter(2)];
  x = x(:);
  y = [y0(1); w(:); y0(2)];
endfunction
