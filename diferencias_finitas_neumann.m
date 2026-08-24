function [x,y] = diferencias_finitas_neumann(f, inter, ya, dyb, N)
  % diferencias finitas
  % f = [p(x); q(x); r(x)]
  % inter = [a b]
  % y0 = y(a)
  % dy0 = y'(b)

  % h = (b-a)/(N+1)
  h = (inter(2) - inter(1))/(N+1);

  % construir sistema
  xi = inter(1) + [1:N+1]*h;
  pqr = f(xi);

  a_i = 2 + h^2 * pqr(2,:);     % diagonal principal
  b_i = -1 + (h/2) * pqr(1,:);  % diagonal derecha
  c_i = -1 - (h/2) * pqr(1,:);  % diagonal izquierda
  d_i = -h^2 * pqr(3,:);

  % condicion dirichlet
  d_i(1) = d_i(1) + (1 + (h/2) * pqr(1,1))*ya;

  % condicion neumann
  d_i(N+1) = d_i(N+1) + (2*h - h^2 * pqr(1,N+1)) * dyb;

  % matriz=
  M = diag(a_i);                         % diagonal principal
  M = M + diag(c_i(2:N+1), -1);   % diagonal izquierda
  M = M + diag(b_i(1:N), 1);  % diagonal derecha

  % seria c_i(N+1)
  M(N+1, N) = -2;

  % resolver sistema
  w = M \ d_i(:);

  % solucion con extremos
  x = [inter(1); xi(:)];
  y = [ya; w(:)];
endfunction
